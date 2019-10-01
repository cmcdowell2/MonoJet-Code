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
    bkg = samples.getSumOfBkg().Clone("Sum of Background")
    data = samples.histo['Data'].Clone("Data")
    

    """texS = TLatex(0.20,0.837173,("#sqrt{s} = 13 TeV, 12.8  fb^{-1}"));
    texS.SetNDC();
    texS.SetTextFont(42);
    texS.SetTextSize(0.040);
    texS.Draw();
    texS1 = TLatex(0.12092,0.907173,"#bf{CMS} : #it{Preliminary}");
    texS1.SetNDC();
    texS1.SetTextFont(42);
    texS1.SetTextSize(0.040);
    texS1.Draw();"""

    ######################################
    xbins = samples.histo['Data'].GetNbinsX(); 
    ybins = samples.histo['Data'].GetNbinsY(); 
    Ratio = samples.histo['Data'].Clone("Ratio Data/MC");
    last_hist = bkg.Clone("Last");
    last = last_hist.Clone("last");
    for xbin in range(1,xbins+1):
        for ybin in range(1, ybins+1):
            stackcontent = last.GetBinContent(xbin,ybin);
            stackerror = last.GetBinError(xbin,ybin);
            datacontent = samples.histo['Data'].GetBinContent(xbin,ybin);
            dataerror = samples.histo['Data'].GetBinError(xbin,ybin);
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
            Ratio.SetBinContent(xbin,ybin,ratiocontent);
            Ratio.SetBinError(xbin,ybin,error);
        


    histos = [bkg, data, Ratio]
    label = {}
    label[bkg.GetName()]={"pre":"MC_","zaxis":[bkg.GetMinimum(),bkg.GetMaximum()]}
    label[data.GetName()]={"pre":"Data_","zaxis":[data.GetMinimum(),data.GetMaximum()]}
    label[Ratio.GetName()]={"pre":"Data-MC_","zaxis":[0.3,1.7]}
    for histo in histos:
        canvas = TCanvas(histo.GetName(), "", 1200, 800)
        canvas.SetMargin(0.15, 0.15, 0.15, 0.08)
        histo.Draw("colz")
        histo.SetStats(0)
        histo.GetZaxis().SetRangeUser(label[histo.GetName()]["zaxis"][0], label[histo.GetName()]["zaxis"][1])
        histo.SetTitle(histo.GetName())
        histo.GetXaxis().SetTitle(samples.name['x'])
        histo.GetYaxis().SetTitle(samples.name['y'])
        histo.GetXaxis().SetTitleOffset(1.7)
        #canvas.Write()
   
        dir = os.getcwd().split("/")[-1]

        file_path = "/afs/hep.wisc.edu/home/mcdowell/public_html/2018Data/" + dir + "/"
        directory = os.path.join(os.path.dirname(file_path), "")

        canvas.SaveAs(directory + label[histo.GetName()]["pre"] + variable + ".pdf");
        canvas.SaveAs(directory + label[histo.GetName()]["pre"] + variable + ".png");


  
