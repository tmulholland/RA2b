#include "TRandom3.h"
#include "TMath.h"


void binomialToy(){

  gStyle->SetOptStat(0);
  //  gStyle->SetTitleSize(1.5,"t");
  //  gStyle->SetTitleX(0.5); //title X location
  //  gStyle->SetTitleY(1.3); //title Y location 


  float njetsRand[13];
  //  float sf[3] = {0.2470, 0.0479, 0.0045};
  float sf[4] = {1, 0.246959,0.0478692,0.00449657};
  
  float njet2, njet3;

  //set probability and gaussian error on probability
  float meanProb = 0.06;
  float deltaProb = 0.06;
  float doLumi = 10.0;
  
  TH1F *histArray[3][2];

  if(doLumi == 10.0){
    histArray[0][0] = new TH1F("histArray","1b extrapolation factor (NJets 7-8)",100,0.35,0.45);
    histArray[1][0] = new TH1F("histArray","2b extrapolation factor (NJets 7-8)",100,0.1,0.17);
    histArray[2][0] = new TH1F("histArray","3b extrapolation factor (NJets 7-8)",100,0.008,0.042);
    
    histArray[0][1] = new TH1F("histArray","1b extrapolation factor (NJets 9+)",100,0.44,0.58);
    histArray[1][1] = new TH1F("histArray","2b extrapolation factor (NJets 9+)",100,0.16,0.29);
    histArray[2][1] = new TH1F("histArray","3b extrapolation factor (NJets 9+)",100,0.02,0.1);
  } else if(doLumi == 3.0){
    histArray[0][0] = new TH1F("histArray","1b extrapolation factor (NJets 7-8)",100,0.3,0.5);
    histArray[1][0] = new TH1F("histArray","2b extrapolation factor (NJets 7-8)",100,0.08,0.2);
    histArray[2][0] = new TH1F("histArray","3b extrapolation factor (NJets 7-8)",100,0.0,0.056);
    
    histArray[0][1] = new TH1F("histArray","1b extrapolation factor (NJets 9+)",100,0.4,0.65);
    histArray[1][1] = new TH1F("histArray","2b extrapolation factor (NJets 9+)",100,0.12,0.35);
    histArray[2][1] = new TH1F("histArray","3b extrapolation factor (NJets 9+)",100,0.0,0.13);
  }
  
  if(doLumi == 10.0) {
    float njets[13] = {0,0,0,0, 666.86, 338.996, 130.993, 41.1673, 10.6336, 2.65273, 0.590898, 0.11701, 0.0246337};
  } else if(doLumi == 3.0) {
    float njets[13] = {0,0,0,0, 199.812, 101.595, 39.3041, 12.3503, 3.19, 0.795817, 0.177269, 0.035103, 0.0073901};
  } else {
    cout << "unknown lumi" << endl;
  }

  
  for(int i=0; i<100000; i++) {

    float prob = getProb(meanProb,deltaProb);
    //    cout << prob << endl;
    
    for(int iter=4; iter<13; iter++){
      njetsRand[iter] = getEvents(njets[iter]);
    }
    
    for(int bJetBin = 1; bJetBin<4; bJetBin++) {
	
      njet2 = sf[bJetBin]*getCorrection(bJetBin,njetsRand,njets,prob,2);
      njet3 = sf[bJetBin]*getCorrection(bJetBin,njetsRand,njets,prob,3);
	
      histArray[bJetBin-1][0]->Fill(njet2);
      //      if((njetsRand[9]>0)&&njetsRand[10]>0)
      histArray[bJetBin-1][1]->Fill(njet3);
      
    }
  }

  float max;
  float cv;
  float error;
  float mean;
  float rms;
  
  TCanvas * c1 = new TCanvas("c1", "c1", 1900, 1000);
  c1->Divide(3,2);
  c1->cd(1);
  histArray[0][0]->Draw();
  histArray[0][0]->SetLineWidth(3);

  max = histArray[0][0]->GetMaximum();
  cv = 0.394919;
  if(doLumi == 10)
    error = 0.0300001;
  if(doLumi == 3)
    error = 0.0548121;
  

  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();
  
  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[0][0]->GetMean();
  rms = histArray[0][0]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

  TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();
  
  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;
  
  //  TCanvas * c2 = new TCanvas("c2", "c2", 600, 400);
  c1->cd(2);
  histArray[1][0]->Draw();
  histArray[1][0]->SetLineWidth(3);

  max = histArray[1][0]->GetMaximum();
  cv = 0.132133;
  if(doLumi == 10)
    error = 0.0208995;
  if(doLumi == 3)
    error = 0.0381918;
  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();

  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[1][0]->GetMean();
  rms = histArray[1][0]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

  TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();

  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;
  
  //  TCanvas * c3 = new TCanvas("c3", "c3", 600, 400);
  c1->cd(3);
  histArray[2][0]->Draw();
  histArray[2][0]->SetLineWidth(3);

  max = histArray[2][0]->GetMaximum();
  cv = 0.0237036;
  if(doLumi == 10)
    error = 0.011977;
  if(doLumi == 3)
    error = 0.0218833;
  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();

  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[2][0]->GetMean();
  rms = histArray[2][0]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

  TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();

  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;
  
  //  TCanvas * c4 = new TCanvas("c4", "c4", 600, 400);
  c1->cd(4);
  histArray[0][1]->Draw();
  histArray[0][1]->SetLineWidth(3);

  max = histArray[0][1]->GetMaximum();
  cv = 0.507533;
  if(doLumi == 10)
    error = 0.0385548;
  if(doLumi == 3)
    error = 0.0704423;
  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();

  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[0][1]->GetMean();
  rms = histArray[0][1]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

  TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();

  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;
  
  //  TCanvas * c5 = new TCanvas("c5", "c5", 600, 400);
  c1->cd(5);
  histArray[1][1]->Draw();
  histArray[1][1]->SetLineWidth(3);

  max = histArray[1][1]->GetMaximum();
  cv = 0.226098;
  if(doLumi == 10)
    error = 0.035762;
  if(doLumi == 3)
    error = 0.0653516;
  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();

  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[1][1]->GetMean();
  rms = histArray[1][1]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

    TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();

  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;

  
  //  TCanvas * c6 = new TCanvas("c6", "c6", 600, 400);
  c1->cd(6);
  histArray[2][1]->Draw();
  histArray[2][1]->SetLineWidth(3);

  max = histArray[2][1]->GetMaximum();
  cv = 0.0565875;
  if(doLumi == 10)
    error = 0.0285926;
  if(doLumi == 3)
    error = 0.052242;
  lowCutLine  = new TLine(cv,0,cv,max);
  lowCutLine->SetLineWidth(3);
  lowCutLine->SetLineStyle(2);
  lowCutLine->Draw();

  TBox *b = new TBox(cv-error,0,cv+error,max);
  b->SetFillColor(1); b->SetFillStyle(3004);
  b->Draw();

  mean = histArray[2][1]->GetMean();
  rms = histArray[2][1]->GetRMS();

  lowCutLine  = new TLine(mean,0,mean,max);
  lowCutLine->SetLineWidth(3);
  //  lowCutLine->SetLineStyle(2);
  lowCutLine->SetLineColor(kRed);
  lowCutLine->Draw();

  TBox *b = new TBox(mean-rms,0,mean+rms,max);
  b->SetFillColor(kRed); b->SetFillStyle(3001);
  b->Draw();

  cout << "toy uncertianty " << 100*rms/mean << "%" << endl;
  cout << "stat uncertianty " << 100*error/cv << "%" << endl;
  
}


// returns random poisson distributed number of events
float getEvents(float events) {
  return (gRandom->Poisson(events));
}

// returns random gaussian distributed number of events
// if pRand < 0, set pRand = 0
float getProb(float p, float dp) {
  float pRand = gRandom->Gaus(p,dp);
  return pRand;
}


float getCorrection(int b, float NJetsRand[], float NJets[], float p, int bin){

  // shouldn't ever have 0 events with 4 jets
  // so there's no need for telling the
  // computer how to divde by 0
  float NJets5o4 = NJetsRand[5]/NJetsRand[4];
  float NJets6o4 = NJetsRand[6]/NJetsRand[4];

  // If no events with 7 jets
  // use mc to get ratios
  if(NJetsRand[7]==0){
    float NJets8o7 = NJets[8]/NJets[7];
  } else {
    float NJets8o7 = NJetsRand[8]/NJetsRand[7];
  }

  // If no events with 9 jets
  // use mc to get ratios
  if(NJetsRand[9]==0){
    float NJets10o9 = NJets[10]/NJets[9];
    float NJets11o9 = NJets[11]/NJets[9];
    float NJets12o9 = NJets[12]/NJets[9];
  } else {
    float NJets10o9 = NJetsRand[10]/NJetsRand[9];
    float NJets11o9 = NJetsRand[11]/NJetsRand[9];
    float NJets12o9 = NJetsRand[12]/NJetsRand[9];
  }
  
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
    corrTerm = 0;
  
  return corrTerm;

}
