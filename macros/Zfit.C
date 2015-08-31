{
  gROOT->Reset();
  gSystem->Load("libRooFit");
  using namespace RooFit ;
  #include <tdrstyle.C>
  //#include "RooKeysPdf.h"

  gROOT->Reset();
  setTDRStyle();
  tdrStyle->SetPadRightMargin(0.05);

  gROOT->ProcessLine(".L setTreeWeight.C");
  gROOT->ProcessLine(".L getCuts.C");    

  bool isData = false;
  
  float doLumi = 382.;
  cout << "chain the samples" << endl;

  TChain dataChain("TreeMaker2/PreSelection") ;

  setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-100to200-RunIISpring15.root",  doLumi, 139.4);
  setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-200to400-RunIISpring15.root",  doLumi, 42.75);
  setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-400to600-RunIISpring15.root",  doLumi, 5.497);
  setTreeWeight("/home/ww/work/data/RA2bTrees/DYJets-600toInf-RunIISpring15.root",  doLumi, 2.21);
  dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-100to200-RunIISpring15.root") ;
  dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-200to400-RunIISpring15.root") ;
  dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-400to600-RunIISpring15.root") ;
  dataChain.Add("/home/ww/work/data/RA2bTrees/DYJets-600toInf-RunIISpring15.root") ;

  int nJetBin = 2;
  int bJetBin = 1;
  int kinBin = -1;

  std::ostringstream ss;
  ss << nJetBin;
  std::string nJ(ss.str());
  std::ostringstream ss;
  ss << bJetBin;
  std::string bJ(ss.str());
  std::ostringstream ss;
  ss << kinBin;
  std::string kin(ss.str());
  
  TString mmpng ="../plots/Zmass/mm_nJ"+nJ+"_bJ"+bJ+"_kin"+kin+".png";
  TString mmpdf ="../plots/Zmass/mm_nJ"+nJ+"_bJ"+bJ+"_kin"+kin+".pdf";
  TString eepng ="../plots/Zmass/ee_nJ"+nJ+"_bJ"+bJ+"_kin"+kin+".png";
  TString eepdf ="../plots/Zmass/ee_nJ"+nJ+"_bJ"+bJ+"_kin"+kin+".pdf";

  TString eepngLoose = "../plots/Zmass/eeLoose.png";
  TString eepdfLoose = "../plots/Zmass/eeLoose.pdf";
  TString mmpngLoose = "../plots/Zmass/mmLoose.png";
  TString mmpdfLoose = "../plots/Zmass/mmLoose.pdf";
  
  TCut looseCutsEE = getCuts(-1,-1,-1,2,false); // baseline to get shapes
  TCut tightCutsEE = getCuts(nJetBin,bJetBin,kinBin,2,false); // bin you measure purity in

  TCut looseCutsMM = getCuts(-1,-1,-1,1,false); // baseline to get shapes
  TCut tightCutsMM = getCuts(nJetBin,bJetBin,kinBin,1,false); // bin you measure purity in

  TH1F *zmassHistoEE = new TH1F("zmassHistoEE","zmassHistoEE", 30, 60, 120);
  TH1F *zmassLooseVLBHistoEE = new TH1F("zmassLooseVLBHistoEE","zmassLooseVLBHistoEE", 30, 60, 120);

  TH1F *zmassHistoMM = new TH1F("zmassHistoMM","zmassHistoMM", 30, 60, 120);
  TH1F *zmassLooseVLBHistoMM = new TH1F("zmassLooseVLBHistoMM","zmassLooseVLBHistoMM", 30, 60, 120);

  dataChain.Project("zmassHistoEE","Zp4.M()",looseCutsEE);
  dataChain.Project("zmassLooseVLBHistoEE","Zp4.M()",tightCutsEE);

  dataChain.Project("zmassHistoMM","Zp4.M()",looseCutsMM);
  dataChain.Project("zmassLooseVLBHistoMM","Zp4.M()",tightCutsMM);

  int maxBinLee = zmassHistoEE->GetMaximumBin();
  float maxLee = zmassHistoEE->GetBinContent(maxBinLee);
  int maxBinLmm = zmassHistoMM->GetMaximumBin();
  float maxLmm = zmassHistoMM->GetBinContent(maxBinLmm);

  x = RooRealVar("x","dimuon mass [GeV]",60,120);
  l = RooArgList(x);

  RooDataHist zmassMM = RooDataHist("zmassMM","zmassMM",l,zmassHistoMM);
  RooDataHist zmassEE = RooDataHist("zmassEE","zmassEE",l,zmassHistoEE);
  RooDataHist zmassLooseVLBmm = RooDataHist("zmassLooseVLBmm","zmassLooseVLBmm",l,zmassLooseVLBHistoMM);
  RooDataHist zmassLooseVLBee = RooDataHist("zmassLooseVLBee","zmassLooseVLBee",l,zmassLooseVLBHistoEE);


  vm = RooRealVar("vm","m in voigtian",91.19,90,92);
  vm.setConstant(kTRUE);
  vg = RooRealVar("vg","g in voigtian",2.64,1.0,3.0);
  vg.setConstant(kTRUE);
  vs = RooRealVar("vs","s in voigtian",1.45,1.0,2.0);
  zfitV = RooVoigtian("voigtian","Voigtian PDF",x,vm,vg,vs);

  c0 = RooRealVar("c0","c0",-0.599,-1,1);
  c1 = RooRealVar("c1","c1",-0.099,-.1,.1);
  c2 = RooRealVar("c2","c2",0.099,-1.,1.);
  bfit = RooChebychev("background","Chebychev PDF",x,RooArgList(c0,c1,c2));
  	
  cbm = RooRealVar("cbm","m in crystal ball",93,88,96);
  cbs = RooRealVar("cbs","s in crystal ball",3.38,2.0,4.0);
  cba = RooRealVar("cba","a in crystal ball",3.38,0.1,5);
  cbn = RooRealVar("cbn","n in crystal ball",10.0,1.0,100.);
  zfitCB = RooCBShape("CBShape","CB PDF",x,cbm,cbs,cba,cbn);
  

  meanCore = RooRealVar("meanCore","meanCore",90.9,80,100);
  widthCore = RooRealVar("widthCore","widthCore",2.0,1.0,5.0);
  meanTail = RooRealVar("meanTail","meanTail",86.5,80,100);
  widthTail = RooRealVar("widthTail","widthTail",5.0,2.0,15.0);
  fracZ = RooRealVar("fracZ","fracZ",0.65,0.40,1.0);
  RooGaussian zfitCore("zfitCore","zfitCore",x,meanCore,widthCore);
  RooGaussian zfitTail("zfitTail","zfitTail",x,meanTail,widthTail);
  RooAddPdf zfitDG("zfitDG","zfitDG",RooArgList(zfitCore,zfitTail),RooArgList(fracZ));

  meanOut = RooRealVar("meanOut","meanOut",90.1,80,100);
  widthOut = RooRealVar("widthOut","widthOut",15.0,3.0,40.0);
  RooGaussian zfitOut("zfitOut","zfitOut",x,meanOut,widthOut);
  fracZout = RooRealVar("fracZout","fracZout",0.05,0.00,0.30);
  RooAddPdf zfitTG("zfitTG","zfitTG",RooArgList(zfitCore,zfitTail,zfitOut),RooArgList(fracZ,fracZout));


  nsig = RooRealVar("nsig","nsig",10000,0,10000000);
  nbkg = RooRealVar("nbkg","nbkg",1000,0,10000000);
  tot = RooRealVar("tot","tot",10000,0,10000000);
  frac = RooRealVar("frac","frac",0.9,0.0,1.0);  

  fsig = RooRealVar("fsig","fsig",0.7,0.2,1.0);


  //
  // Choose fit shape here:
  model = RooAddPdf("model","model",RooArgList(zfitDG,bfit),RooArgList(nsig,nbkg));
  model2 = RooAddPdf("model2","model2",RooArgList(zfitDG,bfit),frac);
  

  //  c0.setConstant(kTRUE);
  //  c1.setConstant(kTRUE);
  //  c2.setConstant(kTRUE);


  // FIRST DO THE MUON PURITY

  // intitialize constants to something reasonable
  c0.setVal(-4.85852e-01);
  c1.setVal(-1.00000e-01);                
  c2.setVal(5.76849e-02);
  fracZ.setVal(5.30594e-01);
  fracZout.setVal(3.00000e-01);
  meanCore.setVal(9.09085e+01);
  meanOut.setVal(8.94648e+01);
  meanTail.setVal(9.04836e+01);
  nbkg.setVal(7.55221e+04);
  nsig.setVal(8.29053e+05);
  widthCore.setVal(1.87054e+00);
  widthOut.setVal(8.03282e+00);
  widthTail.setVal(3.75280e+00);

  model.fitTo(zmassMM);

  float nsigncMM = nsig->getVal();
  float nsigEncMM = nsig->getError();
  float nbkgncMM = nbkg->getVal();
  float nbkgEncMM = nbkg->getError();

  //  TCanvas canvasWfrac("canvasWfrac","canvasWfrac",400,400);
  TCanvas canvas1("canvas1","canvas1",10,10,700,700);
  //  TCanvas canvas("canvas","canvas",10,10,700,700);
  //  canvas.Divide(2,2);
  //  canvas.cd(1);

  int maxBinee = zmassLooseVLBHistoEE->GetMaximumBin();
  float maxee = zmassLooseVLBHistoEE->GetBinContent(maxBinee);
  int maxBinmm = zmassLooseVLBHistoMM->GetMaximumBin();
  float maxmm = zmassLooseVLBHistoMM->GetBinContent(maxBinmm);

  
  lowCutLine  = new TLine(76.19,0.0,76.19,maxee*0.85);
  highCutLine = new TLine(106.19,0.0,106.19,maxee*0.85); 

  lowCutLine2  = new TLine(76.19,0.0,76.19,maxLee*0.85);
  highCutLine2 = new TLine(106.19,0.0,106.19,maxLee*0.85);

  lowCutLine3  = new TLine(76.19,0.0,76.19,maxmm*0.85);
  highCutLine3 = new TLine(106.19,0.0,106.19,maxmm*0.85);

  lowCutLine4  = new TLine(76.19,0.0,76.19,maxLmm*0.85);
  highCutLine4 = new TLine(106.19,0.0,106.19,maxLmm*0.85);


  std::ostringstream ss;
  ss << doLumi;
  std::string L(ss.str());
  
  TString text = "#sqrt{s} = 13TeV, L = "+L+" fb^{-1}";
    
  TPaveText *pt=new TPaveText(0.5,0.92,.98,.98,"brNDC");
  pt->AddText(text);
  pt->SetFillColor(0);

   TPaveText *pt2=new TPaveText(0.01,0.92,.44,.98,"brNDC");


   if(isData)
     pt2->AddText("CMS Preliminary");
   else
     pt2->AddText("CMS Simulation");     

   pt2->SetFillColor(0);


   TPaveText *pt3=new TPaveText(0.18,0.70,.54,.82,"brNDC");
   pt3->AddText("Shapes applied to");
   pt3->AddText("nominal selection");
   //   pt->AddText("PU20bx25");
   pt3->SetFillColor(0);


  xframe = x.frame();
  xframe->GetYaxis()->SetTitle("Events / 2 GeV");
  zmassMM.plotOn(xframe);
  gPad->SetLeftMargin(0.20); //left margin is 15 per cent of the pad width
  gPad->SetTopMargin(0.08); //left margin is 15 per cent of the pad width
  xframe->GetYaxis()->SetTitleOffset(1.7);
  model.plotOn(xframe);
  model.plotOn(xframe,Components("background"),LineStyle(kDashed));

  xframe->SetMaximum(maxLmm*1.20);

  xframe.Draw();
  //  xframe.SaveAs("Zmm_mass.pdf");
   pt->Draw("NB");
   pt2->Draw("NB");


   canvas1.SaveAs(mmpngLoose);
   canvas1.SaveAs(mmpdfLoose);

  lowCutLine4->SetLineWidth(3);
  highCutLine4->SetLineWidth(3);
  lowCutLine4->SetLineStyle(2);
  highCutLine4->SetLineStyle(2);
  //  lowCutLine->SetLineColor(kBlue);
  //  highCutLine->SetLineColor(kBlue);
  lowCutLine4->Draw();
  highCutLine4->Draw();


  c0.setConstant(kTRUE);
  c1.setConstant(kTRUE);
  c2.setConstant(kTRUE);
  meanCore.setConstant(kTRUE);
  widthCore.setConstant(kTRUE);
  meanTail.setConstant(kTRUE);
  widthTail.setConstant(kTRUE);
  meanOut.setConstant(kTRUE);
  widthOut.setConstant(kTRUE);
  fracZ.setConstant(kTRUE);
  fracZout.setConstant(kTRUE);

  x->setRange("all",60,120);
  x->setRange("selection",76.19,106.19);
  //  x->setRange("selection",60.,120.);
  cout << "  integral from 60 to 120 of signal = " << zfitTG.createIntegral(l,"all")->getVal() << endl;
  cout << "  integral from 76.19 to 106.19 of signal = " << zfitTG.createIntegral(l,"selection")->getVal() << endl;
 
  float allMM = zfitTG.createIntegral(l,"all")->getVal();
  float selectionMM = zfitTG.createIntegral(l,"selection")->getVal();
 
  float allBGmm = bfit.createIntegral(l,"all")->getVal();
  float selectionBGmm = bfit.createIntegral(l,"selection")->getVal();

  TCanvas canvas2("canvas2","canvas2",10,10,700,700);
  //  canvas.cd(3);

  model.fitTo(zmassLooseVLBmm);
  xframeLooseVLBmm = x.frame();
  xframeLooseVLBmm->GetYaxis()->SetTitle("Events / 2 GeV");
  zmassLooseVLBmm.plotOn(xframeLooseVLBmm);
  model.plotOn(xframeLooseVLBmm);
  model.plotOn(xframeLooseVLBmm,Components("background"),LineStyle(kDashed));

  xframeLooseVLBmm->SetMaximum(maxmm*1.2);
  xframeLooseVLBmm->SetMinimum(0.);

  gPad->SetTopMargin(0.08); //left margin is 15 per cent of the pad width
  
  xframeLooseVLBmm.Draw();
  lowCutLine3->SetLineWidth(3);
  highCutLine3->SetLineWidth(3);
  lowCutLine3->SetLineStyle(2);
  highCutLine3->SetLineStyle(2);
  //  lowCutLine->SetLineColor(kBlue);
  //  highCutLine->SetLineColor(kBlue);
  lowCutLine3->Draw();
  highCutLine3->Draw();
  pt->Draw("NB");
  pt2->Draw("NB");
  canvas2.SaveAs(mmpng);
  canvas2.SaveAs(mmpdf);

  float nsigLVLBmm = nsig->getVal();
  float nsigELVLBmm = nsig->getError();
  float nbkgLVLBmm = nbkg->getVal();
  float nbkgELVLBmm = nbkg->getError();
  cout << "  Fit finished 7:   nsig = " << nsig << "  nbkg = " << nbkg << endl;
  float allLooseVLBmm = zfitTG.createIntegral(l,Range("all"))->getVal();
  float selectionLooseVLBmm = zfitTG.createIntegral(l,Range("selection"))->getVal();
  //  RooAbsReal* selectionLooseVLBmm = zfitTG.createIntegral(x,NormSet(x),Range("selection"));
  
  float allBGLooseVLBmm = bfit.createIntegral(l,Range("all"))->getVal();
  float selectionBGLooseVLBmm = bfit.createIntegral(l,Range("selection"))->getVal();
  //  RooAbsReal* selectionBGLooseVLBmm = bfit.createIntegral(x,NormSet(x),Range("selection"));

  cout << "  Fit finished 7:   nsigSel = " <<  selectionLooseVLBmm << "  nbkgSel = " << selectionBGLooseVLBmm << endl;

  
  model2.fitTo(zmassLooseVLBmm);
  float looseVLBFractionMM = frac->getVal();
  float looseVLBFractionEmm = frac->getError();

  //NOW DO THE ELECTRON PURITY

  //  canvas.cd(2);

  TCanvas canvas3("canvas3","canvas3",10,10,700,700);

  c0.setConstant(kFALSE);
  c1.setConstant(kFALSE);
  c2.setConstant(kFALSE);
  meanCore.setConstant(kFALSE);
  widthCore.setConstant(kFALSE);
  meanTail.setConstant(kFALSE);
  widthTail.setConstant(kFALSE);
  meanOut.setConstant(kFALSE);
  widthOut.setConstant(kFALSE);
  fracZ.setConstant(kFALSE);
  fracZout.setConstant(kFALSE);

  // intitialize constants to something reasonable
  c0.setVal(-4.85852e-01);
  c1.setVal(-1.00000e-01);                
  c2.setVal(5.76849e-02);
  fracZ.setVal(5.30594e-01);
  fracZout.setVal(3.00000e-01);
  meanCore.setVal(9.09085e+01);
  meanOut.setVal(8.94648e+01);
  meanTail.setVal(9.04836e+01);
  nbkg.setVal(7.55221e+04);
  nsig.setVal(8.29053e+05);
  widthCore.setVal(1.87054e+00);
  widthOut.setVal(8.03282e+00);
  widthTail.setVal(3.75280e+00);

  model.fitTo(zmassEE);
  float nsigncEE = nsig->getVal();
  float nsigEncEE = nsig->getError();
  float nbkgncEE = nbkg->getVal();
  float nbkgEncEE = nbkg->getError();

  xframe = x.frame();
  xframe->GetYaxis()->SetTitle(" ");
  xframe->GetYaxis()->SetTitle("Events / 2 GeV");
  //xframe->GetXaxis()->SetTitle("dielectron mass [GeV]");
  xframe->GetXaxis()->SetTitle("dielectron mass [GeV]");
  // TGaxis::SetMaxDigits(3)

  gPad->SetLeftMargin(0.20); //left margin is 15 per cent of the pad width
  gPad->SetTopMargin(0.08); //left margin is 15 per cent of the pad width
  xframe->GetYaxis()->SetTitleOffset(1.7);
  zmassEE.plotOn(xframe);
  model.plotOn(xframe);
  model.plotOn(xframe,Components("background"),LineStyle(kDashed));

  TGaxis *myX = (TGaxis*)xframe->GetXaxis();

  myX->SetMaxDigits(10);

  xframe->SetMaximum(maxLee*1.20);

  xframe.Draw();

   pt->Draw("NB");

   pt2->Draw("NB");
   
  canvas3.SaveAs(eepngLoose);
  canvas3.SaveAs(eepdfLoose);
  
  lowCutLine2->SetLineWidth(3);
  highCutLine2->SetLineWidth(3);
  lowCutLine2->SetLineStyle(2);
  highCutLine2->SetLineStyle(2);
  //  lowCutLine->SetLineColor(kBlue);
  //  highCutLine->SetLineColor(kBlue);
  lowCutLine2->Draw();
  highCutLine2->Draw();


  c0.setConstant(kTRUE);
  c1.setConstant(kTRUE);
  c2.setConstant(kTRUE);
  meanCore.setConstant(kTRUE);
  widthCore.setConstant(kTRUE);
  meanTail.setConstant(kTRUE);
  widthTail.setConstant(kTRUE);
  meanOut.setConstant(kTRUE);
  widthOut.setConstant(kTRUE);
  fracZ.setConstant(kTRUE);
  fracZout.setConstant(kTRUE);
  
  cout << "  integral from 60 to 120 of signal = " << zfitTG.createIntegral(l,"all")->getVal() << endl;
  cout << "  integral from 76.19 to 106.19 of signal = " << zfitTG.createIntegral(l,"selection")->getVal() << endl;
 
  float allEE = zfitTG.createIntegral(l,"all")->getVal();
  float selectionEE = zfitTG.createIntegral(l,"selection")->getVal();
 
  float allBGee = bfit.createIntegral(l,"all")->getVal();
  float selectionBGee = bfit.createIntegral(l,"selection")->getVal();


  //  canvas.cd(4);

  TCanvas canvas4("canvas4","canvas4",10,10,700,700);

  model.fitTo(zmassLooseVLBee);
  xframeLooseVLBee = x.frame();
  xframeLooseVLBee->GetYaxis()->SetTitle("Events / 2 GeV");
  //  xframeLooseVLBee->SetTitle("Z#rightarrow ee Purity");
  xframeLooseVLBee->GetXaxis()->SetTitle("dielectron mass [GeV]");
  //  xframeLooseVLBee->GetXaxis()->SetTitle("dielectron mass GeV");
  zmassLooseVLBee.plotOn(xframeLooseVLBee);
  model.plotOn(xframeLooseVLBee);
  model.plotOn(xframeLooseVLBee,Components("background"),LineStyle(kDashed));

  xframeLooseVLBee->SetMaximum(maxee*1.20);
  xframeLooseVLBee->SetMinimum(0.);

  xframeLooseVLBee.Draw();
 
   pt->Draw("NB");
   pt2->Draw("NB");


  gPad->SetTopMargin(0.08); //left margin is 15 per cent of the pad width
  
  lowCutLine->SetLineWidth(3);
  highCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  highCutLine->SetLineStyle(2);
  lowCutLine->Draw();
  highCutLine->Draw();

  canvas4.SaveAs(eepng);
  canvas4.SaveAs(eepdf);

  float nsigLVLBee = nsig->getVal();
  float nsigELVLBee = nsig->getError();
  float nbkgLVLBee = nbkg->getVal();
  float nbkgELVLBee = nbkg->getError();
  cout << "  Fit finished 7:   nsig = " << nsig << "  nbkg = " << nbkg << endl;
  float allLooseVLBee = zfitTG.createIntegral(l,"all")->getVal();
  float selectionLooseVLBee = zfitTG.createIntegral(l,"selection")->getVal();
  
  float allBGLooseVLBee = bfit.createIntegral(l,"all")->getVal();
  float selectionBGLooseVLBee = bfit.createIntegral(l,"selection")->getVal();
  
  model2.fitTo(zmassLooseVLBee);
  float looseVLBFractionEE = frac->getVal();
  float looseVLBFractionEee = frac->getError();
  
  cout << endl;   

  cout << "Reporting results: " << endl;
  cout << endl;

  cout << "*******************************DIMUON PURITY*******************************\n" << endl;
  cout << "For full mass window from 60-120" << endl;
  cout << "   No Cut has signal fraction = " << nsigncMM/(nsigncMM+nbkgncMM) << " nsig = " << nsigncMM << "+-" << nsigEncMM << "  nbkg = " << nbkgncMM << "+-" << nbkgEncMM << endl;
  //  cout << "   Loose 1 B   has signal fraction = " << nsigL1B/(nsigL1B+nbkgL1B) << " nsig = " << nsigL1B << "+-" << nsigEL1B << "  nbkg = " << nbkgL1B << "+-" << nbkgEL1B << endl;
  cout << "   LooseVLB   has signal fraction = " << nsigLVLBmm/(nsigLVLBmm+nbkgLVLBmm) << " nsig = " << nsigLVLBmm << "+-" << nsigELVLBmm << "  nbkg = " << nbkgLVLBmm << "+-" << nbkgELVLBmm << endl;
  //cout << "   LooseNoB   has signal fraction = " << nsigLNoB/(nsigLNoB+nbkgLNoB) << " nsig = " << nsigLNoB << "+-" << nsigELNoB << "  nbkg = " << nbkgLNoB << "+-" << nbkgELNoB << endl;
  cout << endl;
  cout << "For mass window of +-15 GeV only" << endl;
  cout << "   mass range of +-15 GeV contains " << selectionMM/allMM << " of all Zmumu events in data " << endl;

  float integralSigMM = ((selectionMM/allMM)*nsigncMM);
  float integralBgMM = ((selectionBGmm/allBGmm)*nbkgncMM);

  cout << "My Purity " <<  selectionLooseVLBmm/(selectionLooseVLBmm+selectionBGLooseVLBmm)  << endl;
  cout << "   fraction events within +-15 GeV that are signal = " << integralSigMM/(integralSigMM+integralBgMM) << endl;

// for error on purity, we'll use the naive binomial estimate scaled up by the difference between the true error from the fit and
// the naive estimate from the full fit range.
  float integralSigLooseVLBmm = ((selectionLooseVLBmm/allLooseVLBmm)*nsigLVLBmm);
  float integralBgLooseVLBmm = ((selectionBGLooseVLBmm/allBGLooseVLBmm)*nbkgLVLBmm);
  float integralSigLooseEVLBmm = ((selectionLooseVLBmm/allLooseVLBmm)*nsigELVLBmm);
  int looseVLBentriesMM = zmassLooseVLBHistoMM->Integral();
  float looseVLBFractionENaiveMM = sqrt( (looseVLBFractionMM*(1-looseVLBFractionMM))/(float)looseVLBentriesMM );
  int looseVLBSelentriesMM = zmassLooseVLBHistoMM->Integral(9,23);
  float looseVLBFractionSelMM = integralSigLooseVLBmm/(integralSigLooseVLBmm+integralBgLooseVLBmm);
  float looseVLBFractionESelNaiveMM = sqrt( (looseVLBFractionSelMM*(1-looseVLBFractionSelMM))/(float)looseVLBSelentriesMM );
  float looseVLBFractionESelMM = looseVLBFractionESelNaiveMM*(looseVLBFractionEmm/looseVLBFractionENaiveMM);
  cout << "looseVLBFraction over whole range = " << looseVLBFractionMM << " +- " << looseVLBFractionEmm << endl;
  cout << "   signal purity in +-15 GeV fit range  = " << integralSigLooseVLBmm/(integralSigLooseVLBmm+integralBgLooseVLBmm) <<
    " +- " << looseVLBFractionESelMM << endl;

  cout << endl;
  cout << "*******************************DIELECTRON PURITY***************************\n" << endl;

  cout << "For full mass window from 60-120" << endl;
  cout << "   No Cut has signal fraction = " << nsigncEE/(nsigncEE+nbkgncEE) << " nsig = " << nsigncEE << "+-" << nsigEncEE << "  nbkg = " << nbkgncEE << "+-" << nbkgEncEE << endl;
  //  cout << "   Loose 1 B   has signal fraction = " << nsigL1B/(nsigL1B+nbkgL1B) << " nsig = " << nsigL1B << "+-" << nsigEL1B << "  nbkg = " << nbkgL1B << "+-" << nbkgEL1B << endl;
  cout << "   LooseVLB   has signal fraction = " << nsigLVLBee/(nsigLVLBee+nbkgLVLBee) << " nsig = " << nsigLVLBee << "+-" << nsigELVLBee << "  nbkg = " << nbkgLVLBee << "+-" << nbkgELVLBee << endl;
  //cout << "   LooseNoB   has signal fraction = " << nsigLNoB/(nsigLNoB+nbkgLNoB) << " nsig = " << nsigLNoB << "+-" << nsigELNoB << "  nbkg = " << nbkgLNoB << "+-" << nbkgELNoB << endl;
  cout << endl;
  cout << "For mass window of +-15 GeV only" << endl;
  cout << "   mass range of +-15 GeV contains " << selectionEE/allEE << " of all Zee events in data " << endl;

  float integralSigEE = ((selectionEE/allEE)*nsigncEE);
  float integralBgEE = ((selectionBGee/allBGee)*nbkgncEE);
  cout << "   fraction events within +-15 GeV that are signal = " << integralSigEE/(integralSigEE+integralBgEE) << endl;

// for error on purity, we'll use the naive binomial estimate scaled up by the difference between the true error from the fit and
// the naive estimate from the full fit range.
  float integralSigLooseVLBee = ((selectionLooseVLBee/allLooseVLBee)*nsigLVLBee);
  float integralBgLooseVLBee = ((selectionBGLooseVLBee/allBGLooseVLBee)*nbkgLVLBee);
  float integralSigLooseEVLBee = ((selectionLooseVLBee/allLooseVLBee)*nsigELVLBee);
  int looseVLBentriesEE = zmassLooseVLBHistoEE->Integral();
  float looseVLBFractionENaiveEE = sqrt( (looseVLBFractionEE*(1-looseVLBFractionEE))/(float)looseVLBentriesEE );
  int looseVLBSelentriesEE = zmassLooseVLBHistoEE->Integral(9,23);
  float looseVLBFractionSelEE = integralSigLooseVLBee/(integralSigLooseVLBee+integralBgLooseVLBee);
  float looseVLBFractionESelNaiveEE = sqrt( (looseVLBFractionSelEE*(1-looseVLBFractionSelEE))/(float)looseVLBSelentriesEE );
  float looseVLBFractionESelEE = looseVLBFractionESelNaiveEE*(looseVLBFractionEee/looseVLBFractionENaiveEE);
  cout << "looseVLBFraction over whole range = " << looseVLBFractionEE << " +- " << looseVLBFractionEee << endl;
  cout << "   signal purity in +-15 GeV fit range  = " << integralSigLooseVLBee/(integralSigLooseVLBee+integralBgLooseVLBee) <<
    " +- " << looseVLBFractionESelEE << endl;

  cout << endl;

}

