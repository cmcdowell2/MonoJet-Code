from ROOT import *
from os import path,system,getcwd,listdir
from samplenames import samplenames


signal_Xsec_file="/nfs_scratch/ekoenig4/MonoZprimeJet/CMSSW_8_0_26_patch1/src/ZprimeTools/monoZprime_XS-2016-correctPDF.txt"
# signal_Xsec_file="/nfs_scratch/ekoenig4/MonoZprimeJet/CMSSW_8_0_26_patch1/src/ZprimeTools/monoZprime_XS-2016-defaultPDF.txt"

class datamc(object):

    def __init__(self,command="noCommand",show=1,lumi=59699):

        #Luminosity
        self.lumi=lumi

        self.show = show

        self.name = 'Xaxis Title'
       
        self.fileDir = "./"
    
        #List of Signal Files and Xsec

        self.Mx_Mv = {}
        self.Mx_Mv_Xsec = {}
        self.signal="null"
        
        #List of Region Data Files
        SignalData_FileNames = ["postMETdata_"];
        SingleEleData_FileNames = ["postSingleEle_"];
        SingleMuData_FileNames = ["postSingleMu_"];
        DoubleEleData_FileNames = ["postDoubleEle_"];
        DoubleMuData_FileNames = ["postDoubleMu_"];

        self.Data_FileNames = {"SignalRegion":SignalData_FileNames,"SingleEle":SingleEleData_FileNames,"SingleMu":SingleMuData_FileNames,"DoubleEle":DoubleEleData_FileNames,"DoubleMu":DoubleMuData_FileNames}
        
        #List of Sample Files and Xsec
        WJets_FileNames = ["postWJets_MLM_","postW100to200_","postW200to400_","postW400to600_","postW600to800_","postW800to1200_","postW1200to2500_","postW2500toInf_"];#"postWJets_MLM",
        WJets_Xsec =      [52940.0,         1395.0         ,407.9          ,57.48            ,12.87           ,5.366            ,1.074             ,0.008001];#50690            ,
        
        ZJets_FileNames = ["postZ100to200_","postZ200to400_","postZ400to600_", "postZ600to800_","postZ800to1200_","postZ1200to2500_","postZ2500toInf_"];
        ZJets_Xsec =      [302.8          ,92.59            ,13.18                ,3.257           ,1.49          ,0.3419          ,0.005146];
        
        GJets_FileNames = ["postGJets40to100_","postGJets100to200_","postGJets200to400_","postGJets400to600_","postGJets600toInf_"];
        GJets_Xsec =      [18700.0             ,8640.0               ,2185.0               ,132.5              ,85.31];
        
        DYJets_FileNames = ["postDY70to100_","postDY100to200_","postDY200to400_","postDY400to600_","postDY600to800_","postDY800to1200_","postDY1200to2500_","postDY2500toInf_"];
        DYJets_Xsec =      [0 , 161.1          ,48.66            ,6.968            ,1.743            ,0.8052           ,0.1933            ,.003468];
        
        TTJets_FileNames = ["postTTJetsDiLept_", "postTTJetsSingleLept_"];
        TTJets_Xsec =      [54.23              , 109.6];
        
        DiBoson_FileNames = ["postWWTo2L2Nu_","postWWToLNuQQ_","postWWTo4Q_","postWW_","postWZ_","postZZ_"];
        DiBoson_Xsec =      [11.08           ,45.99           ,47.73        ,75.8     ,27.6     ,12.14];
        
        QCD_FileNames = ["postQCD100to200_","postQCD200to300_","postQCD300to500_","postQCD500to700_","postQCD700to1000_","postQCD1000to1500_","postQCD1500to2000_","postQCD2000toInf_"];
        QCD_Xsec =      [23700000.0        ,1547000.0         ,322600.0          ,29980.0           ,6334.0             ,1088.0              ,99.11              ,20.23];
        
        self.MC_FileNames = {"WJets":WJets_FileNames,"ZJets":ZJets_FileNames,"GJets":GJets_FileNames,"DYJets":DYJets_FileNames,"TTJets":TTJets_FileNames,"DiBoson":DiBoson_FileNames,"QCD":QCD_FileNames};
        self.MC_Xsec =      {"WJets":WJets_Xsec     ,"ZJets":ZJets_Xsec     ,"GJets":GJets_Xsec     ,"DYJets":DYJets_Xsec     ,"TTJets":TTJets_Xsec     ,"DiBoson":DiBoson_Xsec     ,"QCD":QCD_Xsec};
        self.MC_Color =     {"WJets":kRed-4        ,"ZJets":kMagenta-4      ,"GJets":kGray        ,"DYJets":kTeal-3         ,"TTJets":kOrange-2       ,"DiBoson":kCyan-4         ,"QCD":kBlue-4};
        self.MC_Integral = {"WJets":0     ,"ZJets":0     ,"GJets":0    ,"DYJets":0    ,"TTJets":0     ,"DiBoson":0    ,"QCD":0};

        self.SampleList = ["Data","WJets","ZJets","GJets","DYJets","TTJets","DiBoson","QCD"]

        preRegionData = [".output/postMETdata_0_0_0.root",".output/postSingleEle_0_0_0.root",".output/postSingleMu_0_0_0.root",".output/postDoubleEle_0_0_0.root",".output/postDoubleMu_0_0_0.root"]
        postRegionData =["postMETdata.root","postSingleEle.root","postSingleMu.root","postDoubleEle.root","postDoubleMu.root"] 
        RegionName = ["SignalRegion","SingleEle","SingleMu","DoubleEle","DoubleMu"]

        if command != "noCommand":
            self.region=""
            for i in range(len(RegionName)):
                if path.isfile(preRegionData[i]) or path.isfile(postRegionData[i]): 
		    self.region=RegionName[i]
            if self.region=="":
		    print "No Region Data Files Found, Exiting...";
		    exit()

            if self.region == "SignalRegion":
                self.getSignalXsec(signal_Xsec_file)
                if command[1] == "-1":
                    command.pop(1)
                    self.signal = []
                    mxList = self.Mx_Mv.keys();
		    mxList.sort(key=int);
                    for mx in mxList:
                        mvList = self.Mx_Mv[mx].keys();
  			mvList.sort(key=int)
                        for mv in mvList:
                            self.signal.append("Mx"+mx+"_Mv"+mv)
                else: self.signal = [command[1]]; command.pop(1);
            if self.show == 1:
                print "Running in "+self.region+":"
                print "Plotting at",self.lumi,"pb^{-1}"
            self.HaddFiles()
        
    def initiate(self,variable):
        self.GetVariable(variable)
        self.ScaleHistogram()


    def getSumOfBkg(self):
        sumOfBkg = self.histo["WJets"].Clone()
        for mc in self.MC_Color:
            if (mc != "WJets"):
                sumOfBkg.Add(self.histo[mc])
        return sumOfBkg



    def getSignalXsec(self,xsecFile=signal_Xsec_file,scale=1):
        with open(xsecFile) as f:
            text = f.readlines()
            for line in text:
                if "MonoZPrime_V_Mx" in line:
                    line=line.replace("#","")
                    data=line.replace("MonoZPrime_V_","").split()
                    fn="post"+data[0]
                    xsec=float(line.split(":")[1].split()[0])
                    mx=data[0].split("_")[0].replace("Mx","")
                    mv=data[0].split("_")[1].replace("Mv","")
                    if not mx in self.Mx_Mv:
		    	self.Mx_Mv[mx]={};
		    	self.Mx_Mv_Xsec[mx]={}
                    	self.Mx_Mv[mx][mv]=fn; 
		    	self.Mx_Mv_Xsec[mx][mv]=xsec*scale;
        
    def HaddFiles(self):
        AllFiles={}
        for key, value in self.MC_FileNames.items():
	    AllFiles[key]=value[:]
        AllFiles['Data']=self.Data_FileNames[self.region]
        if self.region=="SignalRegion":
            Mx_Value=self.Mx_Mv.keys();Mx_Value.sort(key=int)
            for mx in Mx_Value:
                Mv_Value=self.Mx_Mv[mx].keys();
		Mv_Value.sort(key=int)
                for mv in Mv_Value:
                    AllFiles["Mx"+mx+"_Mv"+mv]=[self.Mx_Mv[mx][mv]]
                    
        #Hadd files together
        for sample in AllFiles:
            for fn in AllFiles[sample]:
                    if not path.isfile(fn[:-1]+".root") and path:
                        nfile = [f for f in listdir(".output/") if fn in f]
                        if len(nfile) != 0:
                            arg = "hadd -f "+fn[:-1]+".root "
                            for f in nfile:arg+=".output/"+f+" "
                            system(arg)
                
        #Hadd data files together
        # nData=str(len(AllFiles['Data'])-1)
        # if not pat.isfile(AllFiles['Data'][0].split("_")[0]+"_final.root"): system('hadd -f '+AllFiles['Data'][0].split("_")[0]+"_final.root "+AllFiles['Data'][0].split("_")[0]+"_{0.."+nData+"}.root")

    def GetVariableName(self,variable):
        for title in samplenames:
            if title in variable:
               self.name = samplenames[title];
            key = variable.split("_")[-2]
            if key == title:
               self.name = samplenames[title];
               break
                    
    def GetVariable(self,variable):
        self.histo = {}
        self.total = {}

        self.GetVariableName(variable)
        print self.name
        rfile=TFile.Open(self.Data_FileNames[self.region][0][:-1]+".root")
        keys = [keylist.GetName() for keylist in gDirectory.GetListOfKeys()]
        if variable in keys:
	    self.histo['Data']=rfile.Get(variable).Clone();
	    self.histo['Data'].SetDirectory(0)
        else:
	    print "Could not find "+variable+" In "+self.Data_FileNames[self.region][0]+".root, Exiting...";
	    exit()

        if self.region == "SignalRegion":
            for signal in self.signal:
                mx = signal.split("_")[0].replace("Mx","")
                mv = signal.split("_")[1].replace("Mv","")
                rfile=TFile.Open(self.Mx_Mv[mx][mv]+".root")
                keys = [keylist.GetName() for keylist in gDirectory.GetListOfKeys()]
                if variable in keys:
		    hs=rfile.Get(variable).Clone();
		    hs.SetDirectory(0)
                else:
		    print "Could not find "+variable+" In "+self.Mx_Mv[mx][mv]+".root, Exiting...";exit()
                self.histo[signal]=hs
                cutflow=rfile.Get("h_cutflow")
                self.total[signal] = cutflow.GetBinContent(1)
            if not 'Signal' in self.SampleList:
		self.SampleList.append('Signal')

        for sample in self.MC_FileNames:
            self.histo[sample]=[]
            self.total[sample]=[]
            for fn in self.MC_FileNames[sample]:
                rfile=TFile.Open(fn[:-1]+".root")
                keys = [keylist.GetName() for keylist in gDirectory.GetListOfKeys()]
                if variable in keys:
		    hs=rfile.Get(variable).Clone();
		    hs.SetDirectory(0)
                else:
		    print "Could not find "+variable+" In "+fn+".root, Exiting...";
		    exit()
                hs.SetName(variable+"_"+fn)
                self.histo[sample].append(hs)
                cutflow=rfile.Get("h_cutflow")
                self.total[sample].append(cutflow.GetBinContent(1))

    def ScaleHistogram(self):
        print "Scaling Histograms to " + str(self.lumi) + " pb^{-1}"
        for sample in self.SampleList:
            if sample == 'Data':
                integral=(self.histo[sample].Integral())
                space=" "*(15-len(sample))
                if self.show == 1:
		    print "integral of "+sample+space+" here:"+"%.6g" % integral
            elif sample == 'Signal':
                for signal in self.signal:
                    mx = signal.split("_")[0].replace("Mx","")
                    mv = signal.split("_")[1].replace("Mv","")
                    #Scaling = (1/TotalEvents)*Luminosity*NNLO-cross-section
                    scale=(1./self.total[signal])*self.lumi*self.Mx_Mv_Xsec[mx][mv] 
                    self.histo[signal].Scale(scale)
                    integral=(self.histo[signal].Integral())
                    space=" "*(15-len(signal))
                    if self.show == 1:
			print "integral of "+signal+space+" here:"+"%.6g" % integral
            else:
                rawevents = 0
                for i in range(len(self.histo[sample])):
                    if self.MC_FileNames[sample] == "null":
			continue
                    #Scaling = (1/TotalEvents)*Luminosity*NNLO-cross-section
                    rawevents += self.histo[sample][i].Integral()
                    # print self.MC_FileNames[sample][i],self.total[sample][i],self.MC_Xsec[sample][i]
                    scale=(1./self.total[sample][i])*self.lumi*self.MC_Xsec[sample][i]
                    self.histo[sample][i].Scale(scale)
                for i in range(1,len(self.histo[sample])): 
		    self.histo[sample][0].Add(self.histo[sample][i])
                self.histo[sample]=self.histo[sample][0]
                self.histo[sample].SetName(self.histo[sample].GetName().replace(self.MC_FileNames[sample][0],sample))
                integral=(self.histo[sample].Integral())
                self.MC_Integral[sample]=integral
                space=" "*(15-len(sample))
                if self.show == 1:
                    # print "integral of raw"+sample+space+" here:"+"%.6g" % rawevents
                    print "integral of "+sample+space+" here:"+"%.6g" % integral
