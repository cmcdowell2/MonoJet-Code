#!/usr/bin/python

from ROOT import * 
from sys import argv
from sys import path
import DataPlot as plot
import NEWPlot as NEWplot
import os

gROOT.SetBatch(1)

samples = plot.datamc(argv)
samples2 = NEWplot.datamc(argv)

for variable in argv[1:]:
    samples.initiate(variable)
    samples2.initiate(variable)

    c = TCanvas("c", "canvas", 800, 800);
    gStyle.SetOptStat(0);
    gStyle.SetLegendBorderSize(0);
    
    #c.SetLeftMargin(0.15);
    #c.SetLogy();
    #c.cd();
    
    pad1 = TPad("pad1","pad1",0.01,0.25,0.99,0.99);
    pad1.Draw(); pad1.cd();
    pad1.SetLogy();
    pad1.SetFillColor(0); pad1.SetFrameBorderMode(0); pad1.SetBorderMode(0);
    pad1.SetBottomMargin(0.);


    samples.histo['Data'].SetLineWidth(2)
    samples.histo['Data'].SetLineColor(kBlack);
    samples.histo['Data'].SetTitle("");
    samples.histo['Data'].GetXaxis().SetTitle("");
    #samples.histo['Data'].GetXaxis().SetTickLength(0);
    samples.histo['Data'].GetXaxis().SetTitleOffset(999);
    samples.histo['Data'].GetYaxis().SetTitle("");
    #samples.histo['Data'].GetYaxis().SetTickLength(0);
    samples.histo['Data'].GetYaxis().SetTitleOffset(999);
    samples.histo['Data'].SetLineColor(kBlack);
    samples.histo['Data'].SetMarkerStyle(20);
    samples.histo['Data'].SetMarkerColor(kBlack);
    samples.histo['Data'].SetMarkerSize(0.9); 
    samples.histo['Data'].Scale(1/float(samples.histo['Data'].Integral()));

    samples2.histo['Data'].SetLineWidth(2);
    #samples2.histo['Data'].SetLineColor(kWhite);
    samples2.histo['Data'].SetTitle("");
    samples2.histo['Data'].GetXaxis().SetTitle("");
    samples2.histo['Data'].GetXaxis().SetTickLength(0);
    samples2.histo['Data'].GetXaxis().SetLabelOffset(999);
    samples2.histo['Data'].GetYaxis().SetTitle("");
    samples2.histo['Data'].GetYaxis().SetTickLength(0);
    samples2.histo['Data'].GetYaxis().SetLabelOffset(999);
    samples2.histo['Data'].SetLineColor(kRed);
    samples2.histo['Data'].SetMarkerStyle(20);
    samples2.histo['Data'].SetMarkerColor(kRed);
    samples2.histo['Data'].SetMarkerSize(0.9);
    samples2.histo['Data'].Scale(1/float(samples2.histo['Data'].Integral()));



    samples.histo['Data'].Draw("hist same")
    samples2.histo['Data'].Draw('pex0same') 

    leg = TLegend(0.78, 0.80, 0.86, 0.887173, "");
    leg.AddEntry(samples.histo['Data'], "OLD Data", "lp");
    leg.AddEntry(samples2.histo['Data'], "NEW Data", "lp");
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.025);
    leg.Draw();


    lumi_label = "41.5";
    texS = TLatex(0.25, 0.865173, ("#sqrt{s} = 13 TeV, " + lumi_label + " fb^{-1}"));
    texS.SetNDC();
    texS.SetTextFont(42);
    texS.SetTextSize(0.040);
    texS.Draw();
    texS1 = TLatex(0.12092, 0.907173, "#bf{CMS} : #it{Preliminary}");
    texS1.SetNDC();
    texS1.SetTextFont(42);
    texS1.SetTextSize(0.040);
    texS1.Draw();
    
    c.cd();
    pad2 = TPad("pad2","pad2",0.01,0.01,0.99,0.25);
    pad2.Draw(); pad2.cd();
    pad2.SetFillColor(0); pad2.SetFrameBorderMode(0); pad2.SetBorderMode(0);
    pad2.SetTopMargin(0);
    pad2.SetBottomMargin(0.35);


    nbins = samples.histo['Data'].GetNbinsX();
    Ratio = samples.histo['Data'].Clone("Ratio");
    
    for ibin in range(1, nbins + 1):
        data1content = samples.histo['Data'].GetBinContent(ibin);
        data2content = samples2.histo['Data'].GetBinContent(ibin);
        data1error = samples.histo['Data'].GetBinError(ibin);
        data2error = samples2.histo['Data'].GetBinError(ibin);

        ratiocontent = 0;
        error = 0;

        if (data1content != 0 and data2content != 0):
            ratiocontent = (data2content) / data1content
        if (data1content != 0 and data2content != 0):
            error = ratiocontent * ((data1error / data1content)**2 + (data2error / data2content)**2)**(0.5)
        else: 
            error = 2.07

        Ratio.SetBinContent(ibin, ratiocontent);
        Ratio.SetBinError(ibin, error);

    Ratio.GetYaxis().SetRangeUser(0.0, 2.2);
    Ratio.SetStats(0);
    Ratio.GetYaxis().CenterTitle();
    Ratio.SetMarkerStyle(20);
    Ratio.SetMarkerSize(0.7);

    line = TLine(samples2.histo['Data'].GetXaxis().GetXmin(), 1., samples2.histo['Data'].GetXaxis().GetXmax(), 1.);
    line.SetLineStyle(8);
       

    Ratio.Draw("pex0");
    Ratio.GetYaxis().SetLabelSize(0.14);
    Ratio.GetYaxis().SetTitleSize(0.12);
    Ratio.GetYaxis().SetLabelFont(42);
    Ratio.GetYaxis().SetTitleFont(42);
    Ratio.GetYaxis().SetTitleOffset(0.25);
    Ratio.GetYaxis().SetNdivisions(100);
    Ratio.GetYaxis().SetTickLength(0.05);
    
    Ratio.GetXaxis().SetLabelSize(0);
    if (samples.name == "Cutflow"): 
        Ratio.GetXaxis().SetLabelSize(0.14);
    Ratio.GetXaxis().SetTitleSize(0.12);
    Ratio.GetXaxis().SetLabelFont(42);
    Ratio.GetXaxis().SetTitleFont(42);
    Ratio.GetXaxis().SetTitleOffset(0.9);
    Ratio.GetXaxis().SetTickLength(0.05);
    line.SetLineColor(kBlack);
    line.Draw("same");
    
    c.Update();
    samples.histo['Data'].GetYaxis().SetTitle("Events");
    samples.histo['Data'].GetYaxis().SetTitleOffset(1.5);

    nbins = samples.histo['Data'].GetNbinsX();
    xmin = samples2.histo['Data'].GetXaxis().GetXmin();
    xmax = samples2.histo['Data'].GetXaxis().GetXmax();
    xwmin = xmin;
    xwmax = xmax;
   
     
    xaxis = TGaxis(xmin,0,xmax,0,xwmin,xwmax,510);
    xaxis.SetTitle(samples.name);
    if (samples.name == "Cutflow"):
        xaxis.SetTitle("");
    xaxis.SetLabelFont(42);
    xaxis.SetLabelSize(0.10);
    xaxis.SetTitleFont(42);
    xaxis.SetTitleSize(0.12);
    xaxis.SetTitleOffset(1.2);
    xaxis.Draw("SAME");

    if (samples.name == "Cutflow"):
        xaxis.SetTitle("");
        xaxis.SetLabelOffset(999);


    yaxis = TGaxis(xmin,0,xmin,2.2,0,2.2,6,"");
    yaxis.SetTitle("NewData/OldData");
    yaxis.SetLabelFont(42);
    yaxis.SetLabelSize(0.10);
    yaxis.SetTitleFont(42);
    yaxis.SetTitleSize(0.10);
    yaxis.SetTitleOffset(0.35);
    yaxis.Draw("SAME");




    dir = os.getcwd().split("/")[-1]

    file_path = "/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/" + dir + "/"
    directory = os.path.join(os.path.dirname(file_path), "")

    c.SaveAs(directory + "DataComp_" + variable + ".pdf");
    c.SaveAs(directory + "DataComp_" + variable + ".png");
    
