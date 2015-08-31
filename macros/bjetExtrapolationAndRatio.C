/*************************************************
  Plot histograms together, then plot
  their ratio on a histogram below.

  Must provied at least 1 histogram 
  to be the "denomintator", and
  up to 3 histograms to be "numerators".

  All the plots will be shown together in a
  big plot, and below that
  will be a plot showing the ratio of all
  the "numerators" to the "denominator".

  Original Author: 
  Michael Anderson
  March 18, 2009
*************************************************/

#include "TH1F.h"
#include <vector>
#include <cmath>

TCanvas* bjetExtrapolationAndRatio(TH1F* numeratorHist, TH1F* denominatorHist, TString title="", TString xTitle="", TString yTitle="", bool topPlotLogY=false) {  
  vector<TH1F*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist );
  TCanvas* c1 = bjetExtrapolationAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle, topPlotLogY);
  return c1;
}

TCanvas* bjetExtrapolationAndRatio(TH1F* numeratorHist1, TH1F* numeratorHist2, TH1F* denominatorHist, TString title="", TString xTitle="", TString yTitle="", bool topPlotLogY=false) {
  vector<TH1F*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  TCanvas* c1 = bjetExtrapolationAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle, topPlotLogY);
  return c1;
}

TCanvas* bjetExtrapolationAndRatio(TH1F* numeratorHist1, TH1F* numeratorHist2, TH1F* numeratorHist3, TH1F* denominatorHist, TString title="", TString xTitle="", TString yTitle="", bool topPlotLogY=false) {
  vector<TH1F*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  numeratorHistograms.push_back( numeratorHist3 );
  TCanvas *c1 = bjetExtrapolationAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle, topPlotLogY);
  return c1;
}

TCanvas* bjetExtrapolationAndRatio(vector<TH1F*> numeratorHistograms, TH1F* denominatorHist, TString title="", TString xTitle="", TString yTitle="", bool topPlotLogY=false) {

  int numberOfNumeratorHists = numeratorHistograms.size();
  if (numberOfNumeratorHists>3) {
    cout << "Too many histograms for numerator (currently only supports up to 3)" << endl;
    exit;
  }
  if (!denominatorHist) {
    cout << "denominatorHist provided does not exist" << endl;
    exit;
  }

  //*************************************************
  // Variables
  //  bool topPlotLogY = 1;      // 0 = no log; 1= log // now set in function
  TString yTitle2 = "ratio"; // bottom plot y axis title

  vector<int> histColors; 

  histColors.push_back(kRed);  // change colors as you like
  histColors.push_back(kBlack);
  histColors.push_back(kGreen-1);

  vector<int> histSymbols; 
  histSymbols.push_back(kFullSquare);  // change symbols as you like
  histSymbols.push_back(kFullCircle);
  histSymbols.push_back(kFullSquare);

  int histDenominatorColor = kBlue;
  //  int histDenominatorColor = (kGreen+4);

  float defaultRatioYmin = 0.0;
  float defaultRatioYmax = 2.08;
  // END of Variables
  //*************************************************

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,1200,1000);
  c1->Range(0,0,1,1);

  vector<TH1F*> hists;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    hists.push_back( (TH1F*)numeratorHistograms[i]->Clone() );
  }
  TH1F* denominatorHistogram = (TH1F*)denominatorHist->Clone();
  //  denominatorHistogram->Sumw2(); //Sumw2 already called in main program

  // Create ratio histograms
  vector<TH1F*> hist_over_denomHist;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist.push_back( (TH1F*)numeratorHistograms[i]->Clone() );
    hist_over_denomHist[i]->GetTitle();   
    //    hist_over_denomHist[i]->Sumw2(); //Sumw2 already called in main program
    hist_over_denomHist[i]->Divide(denominatorHistogram);
  }


  //*************************************************
  // Bottom plot
  TPad *c1_1 = new TPad("c1_1", "newpad",0.00,0.0.01,0.99,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.3);
  c1_1->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);


  hist_over_denomHist[0]->Draw("e0p X0");
  hist_over_denomHist[0]->SetLineWidth(1);
  hist_over_denomHist[0]->SetLineColor(histColors[0]);
  hist_over_denomHist[0]->SetMarkerStyle(histSymbols[0]);
  hist_over_denomHist[0]->SetMarkerColor(histColors[0]);
  hist_over_denomHist[0]->SetMarkerSize(1.15);
  hist_over_denomHist[0]->SetMinimum(defaultRatioYmin);
  hist_over_denomHist[0]->SetMaximum(defaultRatioYmax);
  hist_over_denomHist[0]->GetYaxis()->SetNdivisions(5);
  hist_over_denomHist[0]->SetTitle(";"+xTitle+";"+yTitle2);
  hist_over_denomHist[0]->GetXaxis()->SetTitleSize(0.16);
  hist_over_denomHist[0]->GetXaxis()->SetLabelSize(0.14);
  hist_over_denomHist[0]->GetXaxis()->SetTitleOffset(0.70);
  hist_over_denomHist[0]->GetYaxis()->SetLabelSize(0.11);
  hist_over_denomHist[0]->GetYaxis()->SetTitleSize(0.14);
  hist_over_denomHist[0]->GetYaxis()->SetTitleOffset(0.28);
  for (int i=1; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist[i]->SetLineWidth(1);
    hist_over_denomHist[i]->SetLineColor(histColors[i]);
    hist_over_denomHist[i]->SetMarkerStyle(histSymbols[i]);
    hist_over_denomHist[i]->SetMarkerSize(1.15);
    hist_over_denomHist[i]->SetMarkerColor(histColors[i]);
    hist_over_denomHist[i]->Draw("same,e0p X0");
  }
  for(int i=1;i<3;i++){
    lowCutLine  = new TLine(i*4+1,defaultRatioYmin,i*4+1,defaultRatioYmax);
    
    lowCutLine->SetLineWidth(3);
    lowCutLine->SetLineStyle(2);
    lowCutLine->Draw();
  }

  // End bottom plot
  //*************************************************


  //*************************************************
  // Top Plot
  c1->cd();
  c1_2 = new TPad("c1_2", "newpad",0.00,0.33,0.99,0.99);
  c1_2->Draw(); 
  c1_2->cd();
  c1_2->SetTopMargin(0.1);
  c1_2->SetBottomMargin(0.01);
  c1_2->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);

  denominatorHistogram->GetXaxis()->SetTitleSize(0.00);
  denominatorHistogram->GetYaxis()->SetLabelSize(0.07);
  denominatorHistogram->GetYaxis()->SetTitleSize(0.08);
  denominatorHistogram->GetYaxis()->SetTitleOffset(0.60);
  denominatorHistogram->SetTitle(title+";;"+yTitle);

  denominatorHistogram->SetLineWidth(2);
  denominatorHistogram->SetLineColor(histDenominatorColor);
  denominatorHistogram->Draw("HIST");
  //  denominatorHistogram->SetLineColor(4);
  //  denominatorHistogram->SetLineWidth(2);
  //  denominatorHistogram->SetLabelSize(0.0);
  TH1F *errorBandHist = (TH1F*)denominatorHistogram->Clone();
  errorBandHist->SetLineColor(histDenominatorColor);
  errorBandHist->SetFillStyle(3018);
  //  errorBandHist->SetFillStyle(1001);

  errorBandHist->SetFillColor(histDenominatorColor);
  errorBandHist->Draw("e2 same");


  //  denominatorHistogram->SetTitleOffset(0.35);
  //  denominatorHistogram->SetTitleSize(0.50,"t");
  gStyle->SetTitleSize(0.08,"t"); 
  gStyle->SetTitleX(0.5); //title X location
  gStyle->SetTitleY(1.0); //title Y location 

  float MaxBin = denominatorHistogram->GetMaximumBin();
  float max = denominatorHistogram->GetBinContent(MaxBin);

  float MinBin = denominatorHistogram->GetMinimumBin();
  float min = denominatorHistogram->GetBinContent(MinBin);
 

  for (int i=0; i<numberOfNumeratorHists; i++) {
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(histColors[i]);
    hists[i]->SetMarkerStyle(histSymbols[i]);
    hists[i]->SetMarkerSize(1.15);
    hists[i]->SetMarkerColor(histColors[i]);
    hists[i]->Draw("same e0p X0");

    float numMaxBin = hists[i]->GetMaximumBin();
    float numMax =  hists[i]->GetBinContent(numMaxBin);
    if( max < numMax )
      max = numMax;
    
    float numMinBin = hists[i]->GetMinimumBin(); 
    float numMin =  hists[i]->GetBinContent(numMinBin);
    if( min > numMin )
      min = numMin;
  }
 


  if (min == 0)
    min = 0.0001;
 
  denominatorHistogram->SetMaximum(max*1.25);
  denominatorHistogram->SetMinimum(min*0.75);

  gStyle->SetOptStat(0);
  //  gStyle->SetErrorX(0);

   for(int i=1;i<3;i++){
     lowCutLine  = new TLine(i*4+1,0.0,i*4+1,max);
     
     lowCutLine->SetLineWidth(3);
     lowCutLine->SetLineStyle(2);
     lowCutLine->Draw();
   }


  TLegend* leg2 = new TLegend(0.5,0.6,0.88,0.88, "","brNDC");  
  leg2->SetBorderSize(2);
  leg2->SetFillStyle(1001); 
  leg2->SetFillColor(kWhite); 

  leg2->AddEntry(denominatorHistogram, "MC", "l");
  leg2->AddEntry(hists[0], "Extrapolation (using NJet:4-6)", "p");
  leg2->AddEntry(hists[1], "Extrapolation with correction", "p");

  leg2->Draw("NB");

  c1_2->SetLogy(topPlotLogY);
  // End bottom plot
  //*************************************************

  return c1;
}
