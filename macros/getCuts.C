// baseline: nJetBin = -1, bJetBin = -1, kinBin = -1, Zcode = 0
// dilepton baseline: nJetBin = -1, bJetBin = -1, kinBin = -1, Zcode = 3

// 3 njet bins: nJetBin = 1, 2, or 3, (-1 for baseline) otherwise
// njet = i will cut on the ith jet multiplicity
// e.g. njet = 7 returns "NJets==7"

// 4 btag bins: bJetBin = 0, 1, 2, or 3, (-1 for baseline)

// kinematic binning: kinBin = 1, 2, 3, 4, 5, or 6, (-1 for baseline)

// sample:
// 0 = sig
// 1 = dimuon
// 2 = dielectron
// 3 = dilepton

TCut getCuts(int nJetBin, int bJetBin, int kinBin, int Zcode) {

  // By default apply mass cut for dilepton.
  // Only not used when doing Z fits.
  return getCuts(nJetBin, bJetBin, kinBin, Zcode, true);   

}

TCut getCuts(int nJetBin, int bJetBin, int kinBin, int Zcode, bool massCut) {

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
    cuts += "(@Muons.size()==2&&@Electrons.size()==0)";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  case 2: // Z->ee
    cuts += "(@Muons.size()==0&&@Electrons.size()==2)";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  case 3: // Z->ll
    cuts += "(@Muons.size()==2&&@Electrons.size()==0)||(@Muons.size()==0&&@Electrons.size()==2)";
    cuts += "DeltaPhi1clean>0.5&&DeltaPhi2clean>0.5&&DeltaPhi3clean>0.3";
    break;
  default: // picked unknown sample
    cuts += "BTags>=99999";
    cout << "******WARNING PICKED UNKNOWN SAMPLE******" << endl;
    cout << "******Plots will be empty******" << endl;
    break;
  }

  if(massCut && (Zcode == 1 || Zcode == 2 || Zcode == 3) )
    cuts += "Zp4.M()>=76.188&&Zp4.M()<=106.188";
  return cuts; 
}
