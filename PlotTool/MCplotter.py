#!/usr/bin/python

from ROOT import * 
from sys import argv
from sys import path
import Plot as plot
import os

gROOT.SetBatch(1)

samples = plot.datamc(argv)
for variable in argv[1:]:
    samples.initiate(variable)

    c = TCanvas("c", "canvas", 800, 800);
    gStyle.SetOptStat(0);
    gStyle.SetLegendBorderSize(0);
    
    c.SetLeftMargin(0.15);
    c.SetLogy();
    c.cd();
    

    
    for mc in samples.MC_Color:
        samples.histo[mc].SetTitle("");
        samples.histo[mc].GetXaxis().SetTitle(samples.name);
        samples.histo[mc].GetXaxis().SetTickLength(0);
        samples.histo[mc].GetXaxis().SetLabelOffset(999);
        samples.histo[mc].GetYaxis().SetTitle("");
        samples.histo[mc].GetYaxis().SetTickLength(0);
        samples.histo[mc].GetYaxis().SetLabelOffset(999);
        samples.histo[mc].SetFillColor(samples.MC_Color[mc]);

    hs_datamc = THStack("hs_datamc","Data/MC comparison");


    hs_order = {}

    if (samples.name == "Cutflow"):
        for key in samples.SampleList:
            if key != "Data":
                hs_order[str(samples.histo[key].GetBinContent(10))] = key

    else:
        for key in samples.MC_Integral:
            hs_order[str(samples.MC_Integral[key])] = key
    keylist = hs_order.keys()
    keylist.sort(key=float)
    for order in keylist:
        hs_datamc.Add(samples.histo[hs_order[order]])
    hs_datamc.SetTitle(samples.name);
    min=0.1;max=pow(10,2.5);
    hs_datamc.SetMinimum(min);
    hs_datamc.SetMaximum(hs_datamc.GetMaximum()*max);


    hs_datamc.Draw("HIST")


    leg = TLegend(0.62, 0.65, 0.86, 0.887, "");
    #leg.AddEntry(samples.histo['Data'], "Data", "lp");
    leg.AddEntry(samples.histo['WJets'],"W#rightarrowl#nu","f");
    leg.AddEntry(samples.histo['DYJets'],"Z#rightarrow ll","F");
    leg.AddEntry(samples.histo['DiBoson'],"WW/WZ/ZZ","F");
    leg.AddEntry(samples.histo['QCD'], "QCD","F");
    leg.AddEntry(samples.histo['TTJets'], "Top Quark", "F");
    leg.AddEntry(samples.histo['GJets'],"#gamma+jets", "F");
    leg.AddEntry(samples.histo['ZJets'],"Z#rightarrow#nu#nu","F");

    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.025);
    leg.Draw();

    texS = TLatex(0.20, 0.865173, ("#sqrt{s} = 13 TeV, 27.6 fb^{-1}"));
    texS.SetNDC();
    texS.SetTextFont(42);
    texS.SetTextSize(0.040);
    texS.Draw();
    texS1 = TLatex(0.12092, 0.907173, "#bf{CMS} : #it{Preliminary}");
    texS1.SetNDC();
    texS1.SetTextFont(42);
    texS1.SetTextSize(0.040);
    texS1.Draw();
   

    dir = os.getcwd().split("/")[-1]

    file_path = "/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/" + dir + "/"
    directory = os.path.join(os.path.dirname(file_path), "")

    c.SaveAs(directory + "MC" + variable + ".pdf");
    c.SaveAs(directory + "MC" + variable + ".png");





    #c.SaveAs((str(variable)+str(".pdf")));
    #c.SaveAs((str(variable)+str(".png")));
    #system((str("mv ")+str(variable)+str(".pdf ")+str("/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/")+str(variable)+str(".pdf")));
    #system((str("mv ")+str(variable)+str(".png ")+str("/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/")+str(variable)+str(".png")));
