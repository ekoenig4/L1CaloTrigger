from ROOT import *
from sys import argv

gROOT.SetBatch(1)

class Jet:
    def __init__(self,iphi,ieta,et):
        self.iphi = iphi
        self.ieta = ieta
        self.et = et
def getJetBox(iphi,ieta,geometry=(7,7),color=kRed):
    dphi,deta=[d/2 for d in geometry]
    box = TBox(iphi-dphi-1,ieta-deta-1,iphi+dphi,ieta+deta)
    box.SetFillStyle(0)
    box.SetLineColor(color)
    box.SetLineWidth(2)
    return box
    return seedlist
def getGCT(ngct=0,grid=True):
    gctmap = {
        '1':( (0,-1),(24,35) ),
        '2':( (24,-1),(48,35) ),
        '3':( (48,-1),(72,35) )
    }
    def GCTBox(geometry):
        (x1,y1),(x2,y2)=geometry
        box = TBox(x1,y1,x2,y2)
        box.SetFillStyle(0)
        box.SetLineColor(kBlack)
        box.SetLineStyle(2)
        box.SetLineWidth(4)
        return box
    def GCTGrid(geometry):
        (x1,y1),(x2,y2)=geometry
        gridlist = []
        for iphi in range(2,25,4):
            for ieta in range(3,36,4):
                rphi,reta = iphi+x1+0.5,ieta+y1-0.5
                gridlist.append( getJetBox(rphi,reta,geometry=(3.0,3.0),color=kBlack) )
        return gridlist
    gctlist,gridlist=[],[]
    if ngct==0:
        gctlist = [ GCTBox( gctmap[str(n)] ) for n in range(1,4) ]
        if grid:
            gridlist = [ GCTGrid( gctmap[str(n)] ) for n in range(1,4) ]
    else:
        gctlist = [ GCTBox( gctmap[str(ngct)] ) ]
        if grid:
            gridlist = [ GCTGrid( gctmap[str(ngct)] ) ]
    return gctlist,gridlist

def plotJets(event,grid=False):
    c=TCanvas('%s_%i' % (event.fname,event.nevent),'%s_%i' % (event.fname,event.nevent),1200,800)
    gStyle.SetOptStat(0);
    gStyle.SetLegendBorderSize(0);
    c.SetGrid()
    event.h_calo.Draw('COLZ')
    event.h_tower.SetLineColor(kRed)
    event.h_tower.SetLineWidth(2)
    event.h_tower.Draw("BOX same")
    event.h_seed.SetLineColor(kRed)
    event.h_seed.SetLineWidth(2)
    event.h_seed.Draw('BOX same')
    boxlist = [ getJetBox(seed.iphi,seed.ieta) for seed in event.jetSeed ]
    for box in boxlist: box.Draw('same')

    if grid:
        gctlist,gridlist = getGCT()
        for gct in gctlist: gct.Draw('same')
        for grid in gridlist:
            for box in grid: box.Draw('same')
        c.gct = (gctlist,gridlist)
    c.boxlist = boxlist
    return c

def compareFiles(e1,e2):
    nj1,nj2 = e1.nJet,e2.nJet
    print '            nj1: %.3f nj2: %.3f ratio: %.3f' % (nj1,nj2,float(nj1)/nj2)
    print 'total jet   et1: %.3f et2: %.3f ratio: %.3f' % (e1.h_seed.Integral(),e2.h_seed.Integral(),e1.h_seed.Integral()/e2.h_seed.Integral())
    print 'total tower et1: %.3f et2: %.3f ratio: %.3f' % (e1.h_tower.Integral(),e2.h_tower.Integral(),e1.h_tower.Integral()/e2.h_tower.Integral())
    print 'total calo  et1: %.3f et2: %.3f ratio: %.3f' % (e1.h_calo.Integral(),e2.h_calo.Integral(),e1.h_calo.Integral()/e2.h_calo.Integral())
    def compare(h1,h2,name):
        ratio=h1.Clone()
        ratio.Divide(h2)
        c=TCanvas(name,name,1200,800)
        gStyle.SetOptStat(0);
        gStyle.SetLegendBorderSize(0);
        c.SetGrid()
        ratio.Draw('COLZ TEXT')
        c.ratio = ratio
        return c
    def diff(h1,h2,name):
        diff=h1.Clone()
        diff.Add(h2,-1)
        c=TCanvas(name,name,1200,800)
        gStyle.SetOptStat(0);
        gStyle.SetLegendBorderSize(0);
        c.SetGrid()
        diff.Draw('COLZ TEXT')
        c.diff = diff
        return c
    return compare(e1.h_seed,e2.h_seed,'seed_ratio_%i'%e1.nevent),diff(e1.h_tower,e2.h_tower,'tower_diff_%i'%e1.nevent)

class Test:
    def __init__(self,fname):
        self.fname = fname
        self.tfile = TFile.Open(fname)
        self.nevents = sum([ key.GetName() == 'seed' for key in self.tfile.GetListOfKeys() ])
        self.store = []
    def getEvent(self,nevent):
        if nevent <= 0 or nevent > self.nevents: return False
        self.nevent = nevent
        self.h_seed = self.tfile.Get('seed;%i'%nevent)
        self.h_seed.SetTitle('Jet Seeds %i' % nevent)
        self.h_norm = self.h_seed.Clone()
        self.h_norm.Divide(self.h_seed)
        self.h_tower= self.tfile.Get('tower;%i'%nevent)
        self.h_tower.SetTitle('Jet Towers %i' % nevent)
        self.h_calo= self.tfile.Get('calo;%i'%nevent)
        self.h_calo.SetTitle('Calo Towers %i'%nevent)
        def getSeeds(seeds):
            seedlist = []
            for iphi in range(1,seeds.GetNbinsX()+1):
                for ieta in range(1,seeds.GetNbinsY()+1):
                    if seeds.GetBinContent(iphi,ieta) > 0:
                        seedlist.append( Jet(iphi,ieta,seeds.GetBinContent(iphi,ieta)) )
            return seedlist
        self.jetSeed= getSeeds(self.h_seed)
        self.nJet = len(self.jetSeed)
        self.store.append( (self.h_seed,self.h_norm,self.h_tower) )

####################
tests = [ Test(fname) for fname in argv[1:] ]
nevents = min( test.nevents for test in tests )
c=TCanvas()
c.Print('test_display.pdf(')
for i in range(nevents):
    print 'Analyzing Event %i' % (i+1)
    for test in tests:
        test.getEvent(i+1)
        plotJets(test).Print('test_display.pdf')

    
    if len(tests) > 1:
        ratio,diff = compareFiles(tests[0],tests[1])
        # ratio.Print('test_display.pdf')
        # diff.Print('test_display.pdf')
c.Print('test_display.pdf)')

from os import system
system('mv test_display.pdf ~/public_html/Trigger/JetAlgo/CMSSW_Plots/')
