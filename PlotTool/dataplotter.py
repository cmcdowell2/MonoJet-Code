#!/usr/bin/python

from ROOT import * 
from sys import argv
from sys import path
import RegionPlot as plot
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
    
    samples.histo['Data'].SetLineWidth(2)
    samples.histo['Data'].SetLineColor(kWhite);
    samples.histo['Data'].SetTitle("");
    samples.histo['Data'].GetXaxis().SetTitle(samples.name);
    samples.histo['Data'].GetXaxis().SetTitleOffset(1.2);
    samples.histo['Data'].GetYaxis().SetTitle("Events");
    samples.histo['Data'].GetYaxis().SetTitleOffset(1.2);
    samples.histo['Data'].SetLineColor(kBlack);
    samples.histo['Data'].SetMarkerStyle(20);
    samples.histo['Data'].SetMarkerSize(0.9);
    
    samples.histo['Data'].Draw('pex0same') 

    leg = TLegend(0.62, 0.72, 0.86, 0.987173, "");
    leg.AddEntry(samples.histo['Data'], "Data", "lp");
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.025);
    leg.Draw();


    lumi_label = "41.5";
    texS = TLatex(0.20, 0.865173, ("#sqrt{s} = 13 TeV, " + lumi_label + " fb^{-1}"));
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

    c.SaveAs(directory + "RegionA" + variable + ".pdf");
    c.SaveAs(directory + "RegionA" + variable + ".png");





    #c.SaveAs((str(variable)+str(".pdf")));
    #c.SaveAs((str(variable)+str(".png")));
    #system((str("mv ")+str(variable)+str(".pdf ")+str("/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/")+str(variable)+str(".pdf")));
    #system((str("mv ")+str(variable)+str(".png ")+str("/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/")+str(variable)+str(".png")));
