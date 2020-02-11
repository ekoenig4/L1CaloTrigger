from ROOT import *
from config import config
import time

def getKey(phi,eta): return str( 100*phi + eta )
def getPhiEta(key): return int(key)/100,int(key)%100
def irange(lo,hi,step=1): return range(lo,hi+1,step)
def BuildCalo(calo):
    calomap = {}
    for ieta in irange(1,calo.GetNbinsY()):
        for iphi in irange(1,calo.GetNbinsX()):
            key = getKey(iphi,ieta)
            calomap[key] = Tower(iphi,ieta,calo.GetBinContent(iphi,ieta))
    return calomap
def keyiter(keylist):
    for key in keylist: yield getPhiEta(key)
    
def getBox(iphi,ieta,geometry=(7,7),color=kBlack):
    dphi,deta=[d/2 for d in geometry]
    box = TBox(iphi-dphi-1,ieta-deta-1,iphi+dphi,ieta+deta)
    box.SetFillStyle(0)
    box.SetLineColor(color)
    box.SetLineWidth(2)
    return box
    
storemap={}
def DrawGrid(gridmap,name=None,pause=False,box=None):
    if not config.debug: return
    if name == None: name = str(len(debug_store))
    if name=='t22':xmax,ymax=22,22
    elif name=='t19':xmax,ymax=19,19
    else:
        xmax = max( iphi for iphi,ieta in keyiter(gridmap) )
        ymax = max( ieta for iphi,ieta in keyiter(gridmap) )

    if name in storemap:
        storemap[name].grid.Delete()
        if hasattr(storemap[name],'center'):
            storemap[name].center.Delete()
        storemap[name].Delete()
    grid = TH2F("grid_%s"%name,";iPhi;iEta",xmax,0,xmax,ymax,0,ymax)
    for i in irange(1,xmax): grid.GetXaxis().SetBinLabel(i,str(i))
    for i in irange(1,ymax): grid.GetYaxis().SetBinLabel(i,str(i))

    for y in irange(1,ymax):
        for x in irange(1,xmax):
            if getKey(x,y) in gridmap:
                grid.SetBinContent(x,y,gridmap[getKey(x,y)].et)
    c = TCanvas("canvas_%s"%name,"")
    c.SetGrid()
    gStyle.SetOptStat(0)
    grid.Draw("COLZ")
    c.grid=grid
    if type(box) is tuple:
        (cx,cy,wx,wy) = box
        center=getBox(cx,cy,geometry=(wx,wy))
        center.Draw('same')
        c.center=center
    elif xmax == 22 and ymax == 22:
        center=getBox(11.5,11.5,geometry=(3.0,3.0))
        center.Draw('same')
        c.center=center
    elif xmax == 19 and ymax == 19:
        center=getBox(10,10)
        center.Draw('same')
        c.center=center
    storemap[getKey(xmax,ymax)]=c
    if pause:
        if config.play: time.sleep(1)
        else: raw_input("%s Paused."%name)
class Tower:
    def __init__(self,iphi,ieta,et):
        self.iphi = iphi
        self.ieta = ieta
        self.et = et
    def __str__(self):
        return 'iPhi: %i iEta: %i Et: %f' % (self.iphi,self.ieta,self.et)
        
class Jet:
    def __init__(self):
        self.iphi = -1
        self.ieta = -1
        self.et = -1
        self.towers = []
    def seed(self,seed):
        if config.debug: print '---Seeding Jet %s' % str(seed)
        self.iphi = seed.iphi
        self.ieta = seed.ieta
        self.seed_et = seed.et
        self.et = seed.et
        self.towers.append(seed)
    def add(self,tower):
        if config.debug and tower.et > 0: print '----Adding Tower %s' % str(tower)
        self.et += tower.et
        self.towers.append(tower)
    def __str__(self):
        return 'iPhi: %i iEta: %i Et: %f' % (self.iphi,self.ieta,self.et)

class Info:
    def __init__(self,fname):
        self.fname = fname
        self.tfile = TFile.Open(fname)
        self.nevents = sum([ key.GetName() == 'seed' for key in self.tfile.GetListOfKeys() ])
    def getCalo(self,i):
        if i < self.nevents:
            calo = self.tfile.Get('calo;%i'%(i+1))
            return BuildCalo(calo)
    def __len__(self): return self.nevents
