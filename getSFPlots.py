from ROOT import *
from sys import argv
import os

files = argv[1:]

fname = {
    "egammaEffi.txt_EGM2D_updatedAll.root":"2018_ElectronRecoSF",
    "2018_ElectronLoose.root":"2018_ElectronLooseID",
    "2018_ElectronTight.root":"2018_ElectronTightID",
    "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta":"2018_MuonTightISO",
    "NUM_LooseRelIso_DEN_LooseID_pt_abseta":"2018_MuonLooseISO",
    "NUM_TightID_DEN_genTracks_pt_abseta":"2018_MuonTightID",
    "NUM_LooseID_DEN_genTracks_pt_abseta":"2018_MuonLooseID"
}

gROOT.SetBatch(1)
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.3f")

for f in files:
    label = f.split("/")[-1]
    c = TCanvas("c","c",1200,800)
    if label ==  "egammaEffi.txt_EGM2D_updatedAll.root" or label ==  "2018_ElectronTight.root" or label ==  "2018_ElectronLoose.root":
        rfile = TFile.Open(f)
        hs = rfile.Get("EGamma_SF2D")
        c.SetLeftMargin(0.17); c.SetRightMargin(0.17); c.SetBottomMargin(0.12)
        hs.Draw("COLZ TEXT")
        for ext in (".png",".pdf"): c.SaveAs(fname[label]+ext)
    if label ==  "RunABCD_SF_ISO.root":
        rfile = TFile.Open(f)
        for name in ("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta","NUM_LooseRelIso_DEN_LooseID_pt_abseta"):
            hs = rfile.Get(name)
            # c.SetLeftMargin(0.2); c.SetRightMargin(0.2); c.SetBottomMargin(0.15)
            hs.Draw("COLZ TEXT")
            for ext in (".png",".pdf"): c.SaveAs(fname[name]+ext)
    if label ==  "RunBCDEF_SF_ID.root":
        rfile = TFile.Open(f)
        for name in ("NUM_TightID_DEN_genTracks_pt_abseta","NUM_LooseID_DEN_genTracks_pt_abseta"):
            hs = rfile.Get(name)
            # c.SetLeftMargin(0.2); c.SetRightMargin(0.2); c.SetBottomMargin(0.15)
            hs.Draw("COLZ TEXT")
            for ext in (".png",".pdf"): c.SaveAs(fname[name]+ext)
        
        
