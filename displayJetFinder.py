from ROOT import *
from sys import argv

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
def getSeeds(seeds,norm=True):
    seedlist = []
    for iphi in range(1,seeds.GetNbinsX()+1):
        for ieta in range(1,seeds.GetNbinsY()+1):
            if seeds.GetBinContent(iphi,ieta) > 0:
                seedlist.append( Jet(iphi,ieta,seeds.GetBinContent(iphi,ieta)) )
                if norm: seeds.SetBinContent(iphi,ieta,1)
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
    
def plotJets(fname):
    tfile = TFile.Open(fname)
    seeds = tfile.Get("seed")
    towers = tfile.Get("tower")
    
    seedlist = getSeeds(seeds)
    
    c=TCanvas(fname,fname,1200,800)
    gStyle.SetOptStat(0);
    gStyle.SetLegendBorderSize(0);
    c.SetGrid()
    towers.Draw('COLZ')
    seeds.SetLineColor(kRed)
    seeds.SetLineWidth(2)
    seeds.Draw('BOX same')
    boxlist = [ getJetBox(seed.iphi,seed.ieta) for seed in seedlist ]
    for box in boxlist: box.Draw('same')

    # gctlist,gridlist = getGCT()
    # for gct in gctlist: gct.Draw('same')
    # for grid in gridlist:
    #     for box in grid: box.Draw('same')
    # c.gct = (gctlist,gridlist)
        
    c.towers = towers
    c.seeds = seeds
    c.boxlist = boxlist
    c.tfile = tfile
    return c

def compareFiles(f1,f2):
    tf1,tf2=TFile.Open(f1),TFile.Open(f2)
    s1,s2=tf1.Get('seed'),tf2.Get('seed')
    t1,t2=tf1.Get('tower'),tf2.Get('tower')
    nj1,nj2=len(getSeeds(s1,norm=False)),len(getSeeds(s2,norm=False))
    print 'nj1: %i nj2: %i ratio: %f' % (nj1,nj2,float(nj1)/nj2)
    print 'set1: %f set2: %f ratio: %f' % (s1.Integral(),s2.Integral(),s1.Integral()/s2.Integral())
    print 'tet1: %f tet2: %f ratio: %f' % (t1.Integral(),t2.Integral(),t1.Integral()/t2.Integral())
    def compare(h1,h2,name):
        ratio=h1.Clone()
        ratio.Divide(h2)
        c=TCanvas(name,name,1200,800)
        gStyle.SetOptStat(0);
        gStyle.SetLegendBorderSize(0);
        c.SetGrid()
        ratio.Draw('COLZ TEXT')
        c.tfiles = (tf1,tf2)
        c.seeds = (h1,h2)
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
        c.tfiles = (tf1,tf2)
        c.seeds = (h1,h2)
        c.diff = diff
        return c
    return compare(s1,s2,'seed_ratio'),compare(t1,t2,'tower_ratio'),diff(t1,t2,'tower_diff')

clist = [ plotJets(fname) for fname in argv[1:] ]
ratio=compareFiles(argv[1],argv[2])
raw_input()
