from ROOT import *
from sys import argv

# gROOT.SetBatch(1)
gStyle.SetOptStat(0)

varmap = {
    'nJet':TH1F("nJet",";nJets",10,0,50),
    'jetPhi':TH1F("jetPhi",";jetPhi",10,0,73),
    'jetEta':TH1F("jetEta",";jetEta",10,0,34),
    'jetEt':TH1F("jetEt",";jetEt",10,0,500)
}

class Info:
    def __init__(self,fname):
        self.fname = fname
        self.tfile = TFile.Open(fname)
        self.tree = self.tfile.Get('tree')
    def get(self,variable,color):
        hs = varmap[variable].Clone('%s_%s' % (self.fname,variable))
        self.tree.Draw('%s>>%s_%s'%(variable,self.fname,variable))
        hs.SetLineColor(color)
        hs.SetLineWidth(2)
        hs.SetMarkerStyle(20)
        hs.SetMarkerSize(1)
        hs.SetFillStyle(0)
        return hs
def compare(fname_num,fname_den):
    num = Info(fname_num)
    den = Info(fname_den)

    output = TFile("jet_efficiency.root","recreate")
    def efficiency(variable):
        print "Analyzing",variable
        h_num = num.get(variable,kBlack)
        h_den = den.get(variable,kRed)
        c = TCanvas(variable,variable)
        c.SetTicks(0,1)
        ratio = TRatioPlot(h_num,h_den)
        ratio.SetH1DrawOpt("pex0")
        ratio.SetH2DrawOpt("hist")
        ratio.SetGraphDrawOpt("p")
        ratio.Draw()
        hi = ratio.GetUpperRefObject()
        hi.GetYaxis().SetTitle("Events")
        

        lo = ratio.GetLowerRefGraph()
        lo.GetYaxis().SetTitle("Old/New")
        lo.GetYaxis().SetRangeUser(0.75,1.35)

        lo.SetMarkerStyle(20)
        lo.SetMarkerSize(1)
        c.Update()
        raw_input()
        c.Write()
    for variable in varmap: efficiency(variable)
    output.Close()

compare(argv[1],argv[2])
