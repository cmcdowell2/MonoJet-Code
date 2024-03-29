#!/usr/bin/python

from ROOT import *
from sys import argv
from sys import path
import Plot as plot
import os

gROOT.SetBatch(1)

samples=plot.datamc(argv)
for variable in argv[1:]:
    samples.initiate(variable)

    c = TCanvas("c", "canvas",800,800);
    gStyle.SetOptStat(0);
    gStyle.SetLegendBorderSize(0);
    #c.SetLeftMargin(0.15);
    #c.SetLogy();
    #c.cd();
    
    pad1 = TPad("pad1","pad1",0.01,0.25,0.99,0.99);
    pad1.Draw(); 
    pad1.cd();
    pad1.SetLogy();
    pad1.SetFillColor(0); 
    pad1.SetFrameBorderMode(0); 
    pad1.SetBorderMode(0);
    pad1.SetBottomMargin(0.);
    
    samples.histo['Data'].SetLineWidth(2)
    samples.histo['Data'].SetLineColor(kWhite);
    samples.histo['Data'].SetTitle("");
    samples.histo['Data'].GetXaxis().SetTitle("");
    samples.histo['Data'].GetXaxis().SetTickLength(0);
    samples.histo['Data'].GetXaxis().SetLabelOffset(999);
    samples.histo['Data'].GetYaxis().SetTitle("");
    samples.histo['Data'].GetYaxis().SetTickLength(0);
    samples.histo['Data'].GetYaxis().SetLabelOffset(999);
    samples.histo['Data'].SetLineColor(kBlack);
    samples.histo['Data'].SetMarkerStyle(20);
    samples.histo['Data'].SetMarkerSize(1.05);

    for mc in samples.MC_Color:
        samples.histo[mc].SetTitle("");
        samples.histo[mc].GetXaxis().SetTitle("");
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
		hs_order[str(samples.histo[key].GetBinContent(11))] = key

    else:
        for key in samples.MC_Integral:
	    hs_order[str(samples.MC_Integral[key])] = key
    keylist = hs_order.keys()
    keylist.sort(key=float)
    for order in keylist:
	hs_datamc.Add(samples.histo[hs_order[order]])
    hs_datamc.SetTitle("");
    min=0.1;max=pow(10,2.5);
    hs_datamc.SetMinimum(min);
    hs_datamc.SetMaximum(hs_datamc.GetMaximum()*max);


    hs_datamc.Draw("HIST")

    samples.histo['Data'].Draw('pex0same')

    if samples.signal != 'null':samples.histo[samples.signal[0]].Draw("HIST SAME")

    #################################################

    leg = TLegend(0.62,0.60,0.86,0.887173,"");
    leg.AddEntry(samples.histo['Data'],"Data","lp");
    if (samples.signal != 'null'): leg.AddEntry(samples.histo[samples.signal[0]], samples.signal[0])   
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

    texS = TLatex(0.20,0.837173,("#sqrt{s} = 13 TeV, 59.7 fb^{-1}"));
    texS.SetNDC();
    texS.SetTextFont(42);
    texS.SetTextSize(0.040);
    texS.Draw();
    texS1 = TLatex(0.12092,0.907173,"#bf{CMS} : #it{Preliminary} (2018)");
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

    ######################################
    nbins = samples.histo['Data'].GetNbinsX();  
    Ratio = samples.histo['Data'].Clone("Ratio");
    last_hist = hs_datamc.GetStack().Last();
    last = last_hist.Clone("last");
    for ibin in range(1,nbins+1):
        stackcontent = last.GetBinContent(ibin);
        stackerror = last.GetBinError(ibin);
        datacontent = samples.histo['Data'].GetBinContent(ibin);
        dataerror = samples.histo['Data'].GetBinError(ibin);
        # print "bin: "+str(ibin)+"stackcontent: "+str(stackcontent)+" and data content: "+str(datacontent)
        ratiocontent=0;
        if(datacontent!=0 and stackcontent != 0):
	    ratiocontent = ( datacontent) / stackcontent
        error=0;
        if(datacontent!=0 and stackcontent != 0): 
	    error = ratiocontent*((dataerror/datacontent)**2 + (stackerror/stackcontent)**2)**(0.5)
        else: 
	    error = 2.07
        # print "bin: "+str(ibin)+" ratio content: "+str(ratiocontent)+" and error: "+str(error);
        Ratio.SetBinContent(ibin,ratiocontent);
        Ratio.SetBinError(ibin,error);
        
    Ratio.GetYaxis().SetRangeUser(0.3,1.7);
    Ratio.SetStats(0);
    
    Ratio.GetYaxis().CenterTitle();
    Ratio.SetMarkerStyle(20);
    Ratio.SetMarkerSize(0.95);
    
    line = TLine(hs_datamc.GetXaxis().GetXmin(), 1.,hs_datamc.GetXaxis().GetXmax(), 1.);
    line.SetLineStyle(8);

    
    Ratio.Draw("pex0");
    Ratio.GetYaxis().SetLabelSize(0.14);
    Ratio.GetYaxis().SetTitleSize(0.12);
    Ratio.GetYaxis().SetLabelFont(42);
    Ratio.GetYaxis().SetTitleFont(42);
    Ratio.GetYaxis().SetTitleOffset(0.25);
    Ratio.GetYaxis().SetNdivisions(100);
    Ratio.GetYaxis().SetTickLength(0.05);
    
    Ratio.GetXaxis().SetLabelSize(0.15);
    Ratio.GetXaxis().SetTitleSize(0.12);
    Ratio.GetXaxis().SetLabelFont(42);
    Ratio.GetXaxis().SetTitleFont(42);
    Ratio.GetXaxis().SetTitleOffset(0.9);
    Ratio.GetXaxis().SetTickLength(0.05);
    line.SetLineColor(kBlack);
    line.Draw("same");
    
    c.Update();
    hs_datamc.GetYaxis().SetTitle("Events");
    hs_datamc.GetYaxis().SetTitleOffset(1.5);
    ###########################################

    nbins = samples.histo['Data'].GetNbinsX();
    xmin = samples.histo['Data'].GetXaxis().GetXmin();
    xmax = samples.histo['Data'].GetXaxis().GetXmax();
    xwmin = xmin;
    xwmax = xmax;
    
    xaxis = TGaxis(xmin,0.3,xmax,0.3,xwmin,xwmax,510);
    xaxis.SetTitle(samples.name);
    xaxis.SetLabelFont(42);
    xaxis.SetLabelSize(0.10);
    xaxis.SetTitleFont(42);
    xaxis.SetTitleSize(0.12);
    xaxis.SetTitleOffset(1.2);
    xaxis.Draw("SAME");

    line = TLine(xmin, 1.,xmax, 1.);
    line.SetLineStyle(8);
    line.SetLineColor(kBlack);
    line.Draw("SAME");
    
    if (samples.name == "Cutflow"):
        xaxis.SetLabelOffset(-999)
        xaxis.SetTitle("");
        label = []
        for i in range(1,nbins+1):
            label.append(TLatex(i-0.5,0.3-0.2,samples.histo['Data'].GetXaxis().GetBinLabel(i)));
	    label[i-1].SetTextSize(0.065);
	    label[i-1].SetTextAngle(-30.);
	    label[i-1].Draw("SAME");
      

    yaxis = TGaxis(xmin,0.3,xmin,1.7,0.3,1.7,6,"");
    yaxis.SetTitle("Data/MC");
    yaxis.SetLabelFont(42);
    yaxis.SetLabelSize(0.10);
    yaxis.SetTitleFont(42);
    yaxis.SetTitleSize(0.12);
    yaxis.SetTitleOffset(0.35);
    yaxis.Draw("SAME");
    
    dir = os.getcwd().split("/")[-1]

    file_path = "/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/" + dir + "/"
    directory = os.path.join(os.path.dirname(file_path), "")

    c.SaveAs(directory + "Data-MC" + variable + ".pdf");
    c.SaveAs(directory + "Data-MC" + variable + ".png");


  
