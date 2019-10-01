//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
////
////To compile using rootcom to an executable named 'analyze':
////$ ./rootcom ZprimeJetsClass analyze
////
////To run, assuming this is compiled to an executable named 'analyze':
////$ ./analyze /hdfs/store/user/uhussain/Zprime_Ntuples/ /cms/uhussain/MonoZprimeJet/CMSSW_8_0_8/src/LightZPrimeAnalysis/JetAnalyzer/test/output.root -1 10000
////Runs over every event in the folder Zprime_Ntuples, reporting progress every 10000 events
////and storing the resulting histograms in the file output.root.
////
//
#define MonoJets_cxx
#include "MonoJets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <vector>

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{ 
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      cout<<"Please enter a valid value for maxEvents (parameter 3)."<<endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      cout<<"Please enter a valid value for reportEvery (parameter 4)."<<endl;
      return 1;
    }
  //const char* file2 = argv[2];

  MonoJets t(argv[1], argv[2], argv[5]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}


double MonoJets::getSF(int mu_index) {
  double muEta_to_use = fabs(muEta->at(mu_index));
  double muPt_to_use = muPt->at(mu_index);
  if      (muPt_to_use >= 120)  muPt_to_use = 119.9;
  else if (muPt_to_use <= 20)  muPt_to_use = 20.1;
  if      (muEta_to_use >= 2.4) muEta_to_use = 2.39;
  
  double tightMuISO_SF_corr = h_tightMuSF_ISO->GetBinContent(h_tightMuSF_ISO->GetXaxis()->FindBin(muPt_to_use),h_tightMuSF_ISO->GetYaxis()->FindBin(muEta_to_use));
  double tightMuID_SF_corr = h_tightMuSF_ID->GetBinContent(h_tightMuSF_ID->GetXaxis()->FindBin(muPt_to_use),h_tightMuSF_ID->GetYaxis()->FindBin(muEta_to_use));
  
  return tightMuISO_SF_corr*tightMuID_SF_corr;
}

double MonoJets::getKfactor(double bosonPt) {
  double EWK_corrected_weight = 1.0*(ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(bosonPt)));
  double NNLO_weight = 1.0*(NNLOCorrection->GetBinContent(NNLOCorrection->GetXaxis()->FindBin(bosonPt)));
  double kfactor = 1;
  if(EWK_corrected_weight!=0 && NNLO_weight!=0)
    kfactor = (EWK_corrected_weight/NNLO_weight);
  else
    kfactor= sample.type == WJets ? 1.21 : 1.23;
  return kfactor;
}

bool MonoJets::inclusiveCut() {
  if (sample.isInclusive)
    return genHT < 100;
  return true;
}

void MonoJets::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal = 0;   


  Long64_t nentries = fChain->GetEntries();
  cout<<"Coming in: "<<endl;
  cout<<"nentries:"<<nentries<<endl;
  Long64_t nentriesToCheck = nentries;   


  double nTotalEvents,nFilters, nHLT, nCRSelection, nMET200, lepMET_MT160, nNoElectrons, nMETcut,nbtagVeto, nDphiJetMET,nJetSelection;
  nTotalEvents = nFilters = nHLT = nCRSelection = nMET200 = lepMET_MT160 = nNoElectrons = nMETcut = nDphiJetMET = nbtagVeto = nJetSelection = 0;

  vector<int> jetveto;

  float dphimin=-99;
  //Event is rejected if it contains a HighPtMuon faking MET

   if (!sample.isData) {
     TFile *weights = TFile::Open("PU_Central.root");
     PU = (TH1D*)weights->Get("pileup");
    
     TFile *f_muSF_ISO = new TFile("RunABCD_SF_ISO.root");
     TFile *f_muSF_ID = new TFile("RunABCD_SF_ID.root");
     h_tightMuSF_ISO = (TH2F*)f_muSF_ISO->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
     h_tightMuSF_ID = (TH2F*)f_muSF_ID->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");
     if (sample.isW_or_ZJet()) {
       //This is the root file with EWK Corrections
       TFile *file = new TFile("kfactors.root");
       if (sample.type == WJets) {
       ewkCorrection = (TH1D*)file->Get("EWKcorr/W");
       NNLOCorrection = (TH1D*)file->Get("WJets_LO/inv_pt");
       } else {
       ewkCorrection = (TH1D*)file->Get("EWKcorr/Z");
       NNLOCorrection = (TH1D*)file->Get("ZJets_LO/inv_pt");
     }
    }
   }
 
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;
  cout<<"Running over "<<nTotal<<" events."<<endl;
  //  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry+=4) {    
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
 
    jetCand     .clear();
    jetveto     .clear();

    double event_weight = 1.;
    int bosonPID;
    double bosonPt;
    bool found = false;
    
    if (!sample.isData) {
      // For each event we find the bin in the PU histogram that corresponds to puTrue->at(0) and store
      // binContent as event_weight
      int bin = PU->GetXaxis()->FindBin(puTrue->at(0));
      event_weight = PU->GetBinContent(bin);
      genWeight > 0.0 ? event_weight*=genWeight : event_weight =0.0;
      if (sample.isW_or_ZJet()) {
	for (int i = 0; i < nMC; i++) {
	  if ((*mcPID)[i] == sample.PID) {
	    found = true;
	    bosonPID = (*mcPID)[i];
	    bosonPt = (*mcPt)[i];
	  }
        }
      }
    }
    
    float metcut= 0.0;
    
    jetCand   = getJetCand(100,2.5,0.8,0.1);
    //CR Variables
    lepindex = -1;
    nTotalEvents++;
    
    if (metFilters == 0 && inclusiveCut()) { 
      nFilters++;
      fillHistos(0,event_weight);
      if (HLTMet>>7&1 == 1 || HLTMet>>8&1 == 1 || HLTMet>>10&1 == 1 || !sample.isData) {
        nHLT++;
	fillHistos(1,event_weight);
	if (jetCand.size() > 0) {
          nJetSelection++;
  	  if (sample.isW_or_ZJet()) 
            event_weight *= getKfactor(bosonPt);
	  vector<int> mulist = muon_veto_tightID(jetCand[0], 20.0);
          vector<int> looseMu = muon_veto_looseID(jetCand[0], 0, 10.0);
	  if (mulist.size() == 1 && looseMu.size() == 1) {
            nCRSelection++;
            fillHistos(2,event_weight);
            if (!sample.isData) 
              event_weight *= getSF(mulist[0]);
            lepindex = mulist[0];
	    jetveto = JetVetoDecision(jetCand[0],lepindex);

	    TLorentzVector lep_4vec;
	    lep_4vec.SetPtEtaPhiE(muPt->at(mulist[0]),muEta->at(mulist[0]),muPhi->at(mulist[0]),muE->at(mulist[0]));
	    lepton_pt = lep_4vec.Pt();

	    TLorentzVector met_4vec;
	    met_4vec.SetPtEtaPhiE(pfMET,0.,pfMETPhi,pfMET);
	    TLorentzVector leptoMET_4vec = lep_4vec + met_4vec;
	    double leptoMET = fabs(leptoMET_4vec.Pt());
	    double leptoMETphi = leptoMET_4vec.Phi();
	    Recoil = leptoMET;

            metcut = (fabs(pfMET - caloMET))/Recoil;
	    if (leptoMET > 250) {
	      nMET200++;
	      fillHistos(3,event_weight);
              vector<int> elelist = electron_veto_looseID(jetCand[0],lepindex,10.);
	      if (elelist.size() == 0) {
	        nNoElectrons++;
                fillHistos(4,event_weight);
		if (metcut < 0.5) {
		  nMETcut++;
		  fillHistos(5,event_weight);
		  float dPhiLepMet = DeltaPhi(muPhi->at(lepindex),pfMETPhi);
		  float lepMET_MT = sqrt(2*muPt->at(lepindex)*pfMET*(1-TMath::Cos(dPhiLepMet)));
		  h_lepMET_MT->Fill(lepMET_MT);
		  if (lepMET_MT < 160) {
		    lepMET_MT160++;
		    fillHistos(6,event_weight);
		    h_metcut->Fill(metcut);
		    if (btagVeto()) {
		      nbtagVeto++;
		      fillHistos(7,event_weight);
		      double minDPhiJetMET_first4 = TMath::Pi();
		      for (int j = 0; j < jetveto.size(); j++) {
                        double dPhiJetMET = DeltaPhi(jetPhi->at(jetveto[j]), pfMETPhi);
			if (dPhiJetMET < minDPhiJetMET_first4) {
			  if (j < 4)
                            minDPhiJetMET_first4 = dPhiJetMET;
		        }   
		      }
		      h_dphimin->Fill(minDPhiJetMET_first4);
		      if (dPhiJetMETcut(jetveto)) {
		        nDphiJetMET++;
		        fillHistos(8,event_weight);
		      }
	            }   
	          }
                }
	      }
            }
          } 
	}
      } 
    }
    tree->Fill();
    
    if (jentry%reportEvery == 0)
    {
      cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<endl;
    }
  }
  h_cutflow->SetBinContent(1,nTotalEvents); 
  h_cutflow->SetBinContent(2,nFilters);
  h_cutflow->SetBinContent(3,nHLT);
  h_cutflow->SetBinContent(4,nJetSelection);
  h_cutflow->SetBinContent(5,nCRSelection);
  h_cutflow->SetBinContent(6,nMET200);
  h_cutflow->SetBinContent(7,nNoElectrons);
  h_cutflow->SetBinContent(8,lepMET_MT160);
  h_cutflow->SetBinContent(9,nMETcut);
  h_cutflow->SetBinContent(10,nbtagVeto);
  h_cutflow->SetBinContent(11,nDphiJetMET);
}

void MonoJets::BookHistos(const char* file2) {
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("MonoJet","MonoJet");
  fileName->cd();
   
  h_cutflow = new TH1D("h_cutflow","h_cutflow",11,0,11);h_cutflow->Sumw2();
  h_cutflow->GetXaxis()->SetBinLabel(1,"Total Events");
  h_cutflow->GetXaxis()->SetBinLabel(2,"metFilters");
  h_cutflow->GetXaxis()->SetBinLabel(3,"Trigger");
  h_cutflow->GetXaxis()->SetBinLabel(4,"GoodJet");
  h_cutflow->GetXaxis()->SetBinLabel(5,"CRSelection"); 
  h_cutflow->GetXaxis()->SetBinLabel(6,"leptoMetCut");
  h_cutflow->GetXaxis()->SetBinLabel(7,"NoElectrons");
  h_cutflow->GetXaxis()->SetBinLabel(8,"lepMET160");
  h_cutflow->GetXaxis()->SetBinLabel(9,"caloMET cut");
  h_cutflow->GetXaxis()->SetBinLabel(10,"B-JetVeto");
  h_cutflow->GetXaxis()->SetBinLabel(11,"DeltaPhiCut");

  float MtBins[51]={180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.,3000.};
  
  float MetBins[45]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1400.,1800.,2000.,2500.};

  float PtBins[54]={100., 120., 140., 160., 180., 200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};
  
  float LeptonPtBins[25] = {20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,250.,300.,350.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.};

  h_metcut  = new TH1F("h_metcut","h_metcut; |pfMET-caloMET|/pfMET", 50,0,1.2);h_metcut->Sumw2();
  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  h_lepMET_MT = new TH1F("h_lepMET_MT","h_lepMET_MT; transverse mass of the lepton-Emiss system",40,0,400);h_lepMET_MT->Sumw2();
  for(int i=0; i<nHisto; i++){
    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    string histname(ptbins);
    h_metFilters[i] = new TH1F(("metFilters"+histname).c_str(),"metFilters",50,0,3000); h_metFilters[i]->Sumw2();
    h_nJets[i]   = new TH1F(("nJets"+histname).c_str(), "nJets;Number of Jets", 10, 0, 10);h_nJets[i]->Sumw2();
    h_pfMETall[i] =  new TH1F(("pfMETall"+histname).c_str(), "pfMET",50,0,2000);h_pfMETall[i] ->Sumw2(); 
    h_pfMET200[i] = new TH1F(("pfMET200"+histname).c_str(), "pfMET",50,170,1500);h_pfMET200[i] ->Sumw2(); 
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "E_{T}^{miss} (GeV)",44,MetBins);h_pfMET[i] ->Sumw2();
    h_pfMETPhi[i] = new TH1F(("pfMETPhi"+histname).c_str(), "pfMETPhi",50,-4,4);h_pfMETPhi[i]->Sumw2();
    h_j1Pt[i]  = new TH1F(("j1pT"+histname).c_str(), "j1pT;p_{T} of Leading Jet (GeV)", 53,PtBins);h_j1Pt[i]->Sumw2();
    h_j1Eta[i] = new TH1F(("j1Eta"+histname).c_str(), "j1Eta; #eta of Leading Jet", 50, -2.5, 2.5);h_j1Eta[i]->Sumw2();
    h_j1Phi[i] = new TH1F(("j1Phi"+histname).c_str(), "j1Phi; #phi of Leading Jet", 50, -3.0, 3.0);h_j1Phi[i]->Sumw2();
    h_j1etaWidth[i] = new TH1F(("j1etaWidth"+histname).c_str(),"j1etaWidh; #eta width of Leading Jet", 50,0,0.25);h_j1etaWidth[i] ->Sumw2();
    h_j1phiWidth[i] = new TH1F(("j1phiWidth"+histname).c_str(),"j1phiWidth; #phi width of Leading Jet", 50, 0,0.5);h_j1phiWidth[i]->Sumw2();
    h_j1Mt[i]  = new TH1F(("j1Mt"+histname).c_str(), "j1Mt;M_{T} of Leading Jet (GeV)", 50,MtBins);h_j1Mt[i]->Sumw2(); 
    h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(),"nVtx;nVtx",70,0,70);h_nVtx[i]->Sumw2(); 
    //CR Histograms
    h_LeptonPt[i] = new TH1F(("h_LeptonPt"+histname).c_str(),"h_LeptonPt",24,LeptonPtBins);h_LeptonPt[i]->Sumw2();
    h_LeptonEta[i] = new TH1F(("h_LeptonEta"+histname).c_str(),"h_LeptonEta",30,-3.0,3.0);h_LeptonEta[i]->Sumw2();
    h_LeptonPhi[i] = new TH1F(("h_LeptonPhi"+histname).c_str(),"h_LeptonPhi",50,-3.1416,3.1416);h_LeptonPhi[i]->Sumw2();
    h_recoil[i] = new TH1F(("h_recoil"+histname).c_str(), "Recoil (GeV)",44,MetBins);h_recoil[i] ->Sumw2();
  }
}

void MonoJets::fillHistos(int histoNumber,double event_weight) {
  if (sample.isData)
     event_weight = 1;
  h_nVtx[histoNumber]->Fill(nVtx, event_weight);
  h_metFilters[histoNumber]->Fill(metFilters, event_weight);
  h_nJets[histoNumber]->Fill(nJet, event_weight);
  h_pfMETall[histoNumber]->Fill(pfMET, event_weight);
  h_pfMET200[histoNumber]->Fill(pfMET, event_weight);
  h_pfMET[histoNumber]->Fill(pfMET, event_weight);
  h_pfMETPhi[histoNumber]->Fill(pfMETPhi, event_weight);
  if (jetCand.size()>0)
  {
    h_j1Pt[histoNumber]->Fill(jetPt->at(jetCand[0]), event_weight);
    h_j1Eta[histoNumber]->Fill(jetEta->at(jetCand[0]), event_weight);
    h_j1Phi[histoNumber]->Fill(jetPhi->at(jetCand[0]), event_weight);
    h_j1Mt[histoNumber]->Fill(jetMt->at(jetCand[0]), event_weight);
  }
  //CR Histograms
  if (lepindex >= 0)
  { 
    h_LeptonPt[histoNumber]->Fill(muPt->at(lepindex), event_weight);
    h_LeptonEta[histoNumber]->Fill(muEta->at(lepindex), event_weight);
    h_LeptonPhi[histoNumber]->Fill(muPhi->at(lepindex), event_weight);
  }
  if (lepton_pt > 0)
    h_recoil[histoNumber]->Fill(Recoil, event_weight);
}

//Function to calculate regular deltaR separate from jet width variable 'dR'
double MonoJets::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deltaeta = fabs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
float MonoJets::DeltaPhi(float phi1, float phi2) {
  float pi = TMath::Pi();
  float dphi = fabs(phi1-phi2);
  if (dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}

//Method to categorize the events based on #of charged Hadrons/tracks in the pencilJet
//int ZprimeJetsClass::getCategory(
vector<int> MonoJets::getJetCand(double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut) {
  vector<int> tmpCand;
  tmpCand.clear();
  for (int p=0;p<nJet;p++)
  {
    bool kinematic = (*jetPt)[p] > jetPtCut && (*jetNHF)[p] < jetNHFCut && (*jetCHF)[p] > jetCHFCut && fabs((*jetEta)[p])<jetEtaCut;
    bool tightJetID = false;
    if ((*jetID)[p]>>0&1 == 1) 
      tightJetID = true;
    if (kinematic && tightJetID) {
      tmpCand.push_back(p);
    }
  }
  return tmpCand;
}

// Method that determines if the jet passes cleaning cuts
vector<int> MonoJets::JetVetoDecision(int jet_index, int mu_index) {
  bool jetVeto=true;
  vector<int> jetindex;

  for (int i = 0; i < nJet; i++) {
     double deltar_mu = deltaR(jetEta->at(i),jetPhi->at(i),muEta->at(mu_index),muPhi->at(mu_index));
     bool tightJetID = false;
     if ((*jetID)[i]>>0&1 == 1) 
       tightJetID = true;
     if (deltar_mu>0.4 && jetPt->at(i) >30.0 && fabs(jetEta->at(i)) < 2.5 && tightJetID)
       jetindex.push_back(i);
  }
  return jetindex;
}

bool MonoJets::btagVeto() {
  bool btagVeto = true;
  for(int i = 0; i < nJet; i++) {
      if(jetPt->at(i) >20.0 && fabs(jetEta->at(i)) < 2.4 && jetCSV2BJetTags->at(i) > 0.84)
	btagVeto = false;
  }
  return btagVeto;
}

bool MonoJets::dPhiJetMETcut(vector<int> jets) {
  //reject jet if it is found within DeltaPhi(jet,MET) < 0.5 
  bool passes = false;
  
  int njetsMax = jets.size();
  //Only look at first four jets (because that's what monojet analysis do)
  if(njetsMax > 4)
    njetsMax = 4;
  int j=0;
  for(;j< njetsMax; j++) {
    if(DeltaPhi((*jetPhi)[j],pfMETPhi) < 0.5)
      break;
  }

  if(j==njetsMax)
    passes = true;

  return passes;  
}

vector<int> MonoJets::electron_veto_tightID(int jet_index, float elePtCut) {
  vector<int> ele_cands;
  ele_cands.clear();
  for (int i = 0; i < nEle; i++) {
    if (eleIDbit->at(i)>>2&1 == 1) {
      if (fabs(eleEta->at(i)) < 2.5) {
        if (elePt->at(i) > elePtCut) {
          if(deltaR(eleEta->at(i),elePhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.3) {
            ele_cands.push_back(i);
          }
        }
      }
    }
  }         
  return ele_cands;
}

vector<int> MonoJets::muon_veto_tightID(int jet_index, float muPtCut) {
  vector<int> mu_cands;
  mu_cands.clear();

  for (int i = 0; i < nMu; i++) {
    if(muIDbit->at(i)>>3&1 == 1 && muIDbit->at(i)>>9&1==1) {
      if (fabs(muEta->at(i)) < 2.4) {
	if(muPt->at(i) > muPtCut) {
	  if(deltaR(muEta->at(i),muPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.4) {
            mu_cands.push_back(i);
	  }
	}
      }
    }
  }
  return mu_cands;
}

vector<int> MonoJets::electron_veto_looseID(int jet_index, int mu_index, float elePtCut) {
  vector<int> ele_cands;
  ele_cands.clear();

  for(int i = 0; i < nEle; i++) {
    if(eleIDbit->at(i)>>0&1 == 1) {
      if (fabs(eleEta->at(i)) < 2.5) {
	if(elePt->at(i) > elePtCut) {
	  if(deltaR(eleEta->at(i),elePhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.3)
          {
            ele_cands.push_back(i);
	  }
	}
      }
    }
  }
  return ele_cands;
}

//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5
vector<int> MonoJets::muon_veto_looseID(int jet_index, int ele_index, float muPtCut) {
  vector<int> mu_cands;
  mu_cands.clear();

  for(int i = 0; i < nMu; i++) {
    if(muIDbit->at(i)>>0&1==1) {
      if (fabs(muEta->at(i)) < 2.4) {
	if(muPt->at(i) > muPtCut) {
	  if(deltaR(muEta->at(i), muPhi->at(i), jetEta->at(jet_index), jetPhi->at(jet_index)) > 0.4) {
            mu_cands.push_back(i);
	  }
	}
      }
    }
  }
  return mu_cands;
}

