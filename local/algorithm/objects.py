from ROOT import *
from config import config

def getKey(phi,eta): return str( 100*phi + eta )
def getPhiEta(key): return int(key)/100,int(key)%100
def irange(lo,hi,step=1): return range(lo,hi,step)
def BuildCalo(calo):
    calomap = {}
    for iphi in irange(1,calo.GetNbinsX()+1):
        for ieta in irange(1,calo.GetNbinsY()+1):
            calomap[getKey(iphi,ieta)] = Tower(iphi,ieta,calo.GetBinContent(iphi,ieta))
    return calomap

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
        if config.debug: print '--Seeding Jet %s' % str(seed)
        self.iphi = seed.iphi
        self.ieta = seed.ieta
        self.seed_et = seed.et
        self.et = seed.et
        self.towers.append(seed)
    def add(self,tower):
        if config.debug and tower.et > 0: print '---Adding Tower %s' % str(tower)
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
